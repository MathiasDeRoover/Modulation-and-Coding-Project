clear
close all
clc
addpath(genpath(pwd))

%% Initialization
N = 100;
SNR = 20;
bps = 1;

bitStream = CreateBitStream(N,bps);

%% Send bitStream through ideal channel (ignore LDPC for now)
% bitStreamRec = IdealChannel_exec( bitStream,SNR,'BPSK','no_det' );

%% Starting configurations
ftaps   = 80;            % Amount of causal (and non causal) filter taps
symRate = 2*1e6;         % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T       = 1/symRate;     % Symbol period
modulation  = 'psk';

%% Symbol mapping, upsampling and shaping filter
symStream = mapping(bitStream, bps, modulation);

M = 2;                                                 % UpSample factor
fs = symRate * M;                                       % Sample Frequency
supStream = upsample(symStream,M);

beta = 0.3;                                             % Roll off factor
H = RRCFilter( T,fs,beta,ftaps )';                      % Creating the window
h = ifft(H,'symmetric');                                % Transforming to the time domain
h = h/h(1);                                             % Normalizing the window in the time domain
H = fft(h);                                             % Transforming the normalized window back to the frequency domain

G = sqrt(H);                                            % G is the square root of H so that G*G = H (after the convolutions)
g = fftshift(ifft(G,'symmetric'));                      % Transform G to the time domain. Fftshift is needed to get a proper raised cosine

sgStream = conv(supStream,g,'same');                           % Windowing upStream in the frequencyDomain with g(t)

%% Ignore noise for now
% SignalEnergy = (trapz(abs(sgStream).^2))*(1/fs);        % Total signal energy (this is given in slides) %Trapz is integral approx.
% Eb = SignalEnergy / (N*bps);                            % Energy for a single bit? see slides
% Eb = Eb/2;                                              % The power of a bandpass signal is equal to the power of its complex envelope divided by 2
% EbN0 = db2mag(20*2);                                    % Variable SNR; factor 2 because we are working with powers, not amplitudes: powerdB = 10*log10(power) vs amplitudedB = 20*log10(amplitude)
% N0 = Eb/EbN0;
% NoisePower = 2*N0*fs;
% 
% noise = sqrt(NoisePower/2) * (randn(numel(sgStream),1) + 1i*randn(numel(sgStream),1));
% sgStream = sgStream + noise;



%% Receiver
gmin=fliplr(g);                                         % Converting g(t) to g(-t) to get matched filter
% Shift gmin over eps*T to simulate sampling time offset. 
% r(t) = sum(I(n)*g(t-nT-eps*T))

filter_t     = ((-floor(numel(h)/2):1:(floor(numel(h)/2)-1))*1/fs).';
eps          = 0.5;                                 % between [0, 1]
plot_t_shift = filter_t + eps/fs;

% We need to find the values of gmin at the locations corresponding to
% plot_t_shift

gminshift = interp1(filter_t,gmin,plot_t_shift);
gminshift(end) = gminshift(1);                       % Interpolation is not possible for last value

% Show result
figure
plot(filter_t,gmin)
hold on
plot(plot_t_shift,gmin)
plot(plot_t_shift,gminshift,'*')
hold off
legend('original filter','shifted filter','original filter evaluated at shifted points')


switch modulation                                       % Windowing downStream in the frequencyDomain with g(-t)
    case 'pam'
        sggStream = real(conv(sgStream,gminshift));          % Taking real because noise is complex and pam signal is real.
    otherwise
        sggStream = conv(sgStream,gminshift,'same');                % Noise and qam signal are complex.
end


sggStreamCorrect = conv(sgStream,gmin,'same');
figure
stem(real(sggStream(1:10)))
hold on
stem(real(sggStreamCorrect(1:10)))
hold off
legend('Wrong symbols','Correct symbols')
title('Comparison of correct symbols with shifted sampled symbols')

%% 2nd Method: get symbol values at sample points shifted by eps*T
symTime = 0:1/fs:(N*M-1)/fs;
sampleTime = symTime + eps/fs;
sampleValues = interp1(symTime,sggStreamCorrect,sampleTime);

figure
stem(real(sggStream(1:10)))
hold on
stem(real(sampleValues(1:10)))
hold off
legend('Shifted filter symbols','Shifted sample symbols')
title('Comparision of both methods')

%% Implement Gardner algorithm
% y_eps(n) is the matched filter output sampled at nT + eps*T
epsEst = zeros(N,1);
epsEst2 = zeros(N,1);
epsEst3 = zeros(N,1);
K = 0.05;
for n=1:N-1
    epsEst(n+1) = epsEst(n) + 2*K*real(sggStream(2*n)*(conj(sggStream(2*n+1))-conj(sggStream(2*n-1))));
    epsEst2(n+1) = epsEst2(n) + 2*K*real(sampleValues(2*n)*(conj(sampleValues(2*n+1))-conj(sampleValues(2*n-1))));
end
epsEst(end)
epsEst2(end)
% epsEst should converge to the real shift value


rmpath(genpath(pwd))