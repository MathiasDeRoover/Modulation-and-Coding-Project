clear
close all
clc
addpath(genpath(pwd))

%% Initialization
N = 1000;
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
% supStream = symStream;

beta = 0.3;                                             % Roll off factor
H = RRCFilter( T,fs,beta,ftaps )';                      % Creating the window
h = ifft(H,'symmetric');                                % Transforming to the time domain
h = h/h(1);                                             % Normalizing the window in the time domain
H = fft(h);                                             % Transforming the normalized window back to the frequency domain

G = sqrt(H);                                            % G is the square root of H so that G*G = H (after the convolutions)
g = fftshift(ifft(G,'symmetric'));                      % Transform G to the time domain. Fftshift is needed to get a proper raised cosine

sgStream = conv(supStream,g,'same');                           % Windowing upStream in the frequencyDomain with g(t)

%% Ignore noise for now

%% Receiver
gmin=fliplr(g);                                         % Converting g(t) to g(-t) to get matched filter
% Shift gmin over eps*T to simulate sampling time offset. 
% r(t) = sum(I(n)*g(t-nT-eps*T))

plot_t       = ((-floor(numel(h)/2):1:(floor(numel(h)/2)-1))*1/fs).';
eps          = 0.5;                                              % between [0, 1]
plot_t_shift = plot_t + eps/fs;

% We need to find the values of gmin at the locations corresponding to
% plot_t_shift

gminshift = interp1(plot_t,gmin,plot_t_shift);
gminshift(end) = gminshift(1);                       % Interpolation is not possible for last value

% Show result
figure
plot(plot_t,gmin)
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


% sggStreamCorrect = conv(sgStream,gmin,'same');
% figure
% stem(real(sggStream))
% hold on
% stem(real(sggStreamCorrect))
% hold off


% Instead of dropping samples here, convolutions have been updated with the
% 'same' argument. Otherwise the code does not work for M = 2.
% sggStream = sggStream(2*ftaps+1:end-2*ftaps);           % Dropping access data that originates from convolutions

%% Implement Gardner algorithm
% y_eps(n) is the matched filter output sampled at nT + eps*T
epsEst = zeros(N,1);
K = 0.05;
for n=1:N-1
    epsEst(n+1) = epsEst(n) + K*real(sggStream(2*n)*(conj(sggStream(2*n+1))-conj(sggStream(2*n-1)))); 
end
% epsEst(end)
% epsEst should converge to the real shift value





rmpath(genpath(pwd))