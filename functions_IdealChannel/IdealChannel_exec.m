function [ output_stream ] = IdealChannel_exec( in_stream,SNR,modu,option,uncoded)
%IdealChannel_exec A function to execute the ideal channel calculations
%   Takes an input stream, SNR and modulation method.

%% Starting configurations
ftaps   = 80;            % Amount of causal (and non causal) filter taps
symRate = 2*1e6;         % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T       = 1/symRate;     % Symbol period

%% Tranceiver
switch modu
    case 'BPSK'
        modulation  = 'psk';
        bps         = 1;
    case 'QPSK'
        modulation  = 'psk';
        bps         = 2;
    case '16QAM'
        modulation  = 'qam';
        bps         = 4;
    case '64QAM'
        modulation  = 'qam';
        bps         = 6;
    otherwise
        error([modu,' is an unsupported constellation'])
end

bitStream = in_stream;
if mod(numel(bitStream),bps) ~= 0
    error('Number of bits not a multiple of bps');
end
N = numel(bitStream)/bps;

symStream = mapping(bitStream, bps, modulation);

M = 10;                                                 % UpSample factor
fs = symRate * M;                                       % Sample Frequency
supStream = upsample(symStream,M);

beta = 0.3;                                             % Roll off factor
H = RRCFilter( T,fs,beta,ftaps )';                      % Creating the window
h = ifft(H,'symmetric');                                % Transforming to the time domain
h = h/h(1);                                             % Normalizing the window in the time domain
H = fft(h);                                             % Transforming the normalized window back to the frequency domain

G = sqrt(H);                                            % G is the square root of H so that G*G = H (after the convolutions)
g = fftshift(ifft(G,'symmetric'));                      % Transform G to the time domain. Fftshift is needed to get a proper raised cosine

sgStream = conv(supStream,g);                           % Windowing upStream in the frequencyDomain with g(t)

%% Ideal Channel
% SignalEnergy = (trapz(abs(sgStream).^2))*(1/fs);       % Total signal energy (this is given in slides) %Trapz is integral approx.
SignalEnergy = (trapz(abs(uncoded).^2))*(1/fs); 
Eb = SignalEnergy / (N);                                % Energy for a single bit? see slides
Eb = Eb/2;                                              % The power of a bandpass signal is equal to the power of its complex envelope divided by 2
EbN0 = db2mag(SNR*2);                                   % Variable SNR; factor 2 because we are working with powers, not amplitudes: powerdB = 10*log10(power) vs amplitudedB = 20*log10(amplitude)
N0 = Eb/EbN0;
NoisePower = 2*N0*fs;

noise = sqrt(NoisePower/2) * (randn(numel(sgStream),1) + 1i*randn(numel(sgStream),1));
sgStream = sgStream + noise;


%% Receiver
gmin=fliplr(g);                                         % Converting g(t) to g(-t) to get matched filter
switch modulation                                       % Windowing downStream in the frequencyDomain with g(-t)
    case 'pam'
        sggStream = real(conv(sgStream,gmin));          % Taking real because noise is complex and pam signal is real.
    otherwise
        sggStream = conv(sgStream,gmin);                % Noise and qam signal are complex.
end

sggStream = sggStream(2*ftaps+1:end-2*ftaps);           % Dropping access data that originates from convolutions
shsStream = sggStream(1:M:end);                         % Sampling at nT

switch option
    case 'no_det'
        output_stream = shsStream;
    case 'det'
        output_stream = demapping(shsStream, bps, modulation);      % Demapping
    otherwise
        error('Option not defined')

end

