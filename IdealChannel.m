%% Starting configurations
clear
close all
clc

N       = 1000;          % Amount of Symbols
ftaps   = 80;            % Amount of causal (and non causal) filter taps
modu     = 'QPSK';      % Choose type of modulation (BPSK, QPSK, 16QAM, 64QAM)
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

bitStream = CreateBitStream(N,bps);
if mod(numel(bitStream),bps) ~= 0
    error('Number of bits not a multiple of bps');
end

symStream = mapping(bitStream, bps, modulation);

M = 4;                                                  % UpSample factor
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

%%% Plotting the window
plot_f = -fs/2:fs/2/ftaps:fs/2-fs/2/ftaps;                  % Frequency axis window 
plot_t = (-floor(numel(h)/2):1:floor(numel(h)/2)-1)*1/fs;   % Time axis window
% plot_f = -fs/2:fs/2/ftaps:fs/2;                         % Frequency axis window 
% plot_t = (-floor(numel(h)/2):1:floor(numel(h)/2))*1/fs; % Time axis window
plot_H = fftshift(real(H));                             % H ready to plot
plot_h = fftshift(ifft(H,'symmetric'));                 % h ready to plot
plot_g = fftshift(ifft(G,'symmetric'));                 % g ready to plot
plot_gmin = fliplr(plot_g);                             % g(-t) ready to plot

figure()
subplot(2,2,1) %plot1
plot(plot_f,plot_H);
title('H in f domain')
xlabel('Frequency [Hz]')
ylabel('H(f)')
subplot(2,2,2) %plot2
hold on
plot(plot_t,plot_h);
stem(plot_t(1:M:end),plot_h(1:M:end))
title('H in t domain')
xlabel('time [s]')
ylabel('h(t)')
hold off
subplot(2,2,3) %plot3
plot(plot_t,plot_g);
title("g(t)")
xlabel('time [s]')
ylabel('g(t)')
subplot(2,2,4) %plot4
plot(plot_t,plot_gmin);
title("g(-t)")
xlabel('time [s]')
ylabel('g(-t)')
%%%

% %%Design of filter in time domain
% t_plot=plot_t;
% htest = sin(pi*t_plot/T)/pi./t_plot*T.*cos(pi*beta*t_plot./T)./(1-4*beta.^2*t_plot.^2./T^2);
% htest=htest/max(htest);
% figure
% subplot(1,2,1)
% stem(t_plot,htest);
% xlabel('t')
% ylabel('h(t)')
% subplot(1,2,2)
% htest=[0 htest];%To get max index right
% stem(htest(1:M:end));
% xlabel('t')
% ylabel('h(t)')

% figure 
% plot(plot_t,plot_h);
% hold on
% stem(plot_t,plot_h);
% stem(plot_t(1:M:end),plot_h(1:M:end))
% title('H in t domain')
% xlabel('time [s]')
% ylabel('h(t)')


%% Ideal Channel
SignalEnergy = (trapz(abs(sgStream).^2))*(1/fs);        % Total signal energy (this is given in slides) %Trapz is integral approx.
Eb = SignalEnergy / (N*bps);                            % Energy for a single bit? see slides
%Eb = Eb/2;%Why? We work in baseband so...?              % The power of a bandpass signal is equal to the power of its complex envelope divided by 2
% EbN0 = db2mag(40);                                      % Variable SNR ?
EbN0 = db2mag(20*2);
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
recStream = demapping(shsStream, bps, modulation);      % Demapping

%%% Plotting the received signals
figure
hold on
plot(real(shsStream),imag(shsStream),'*')
plot(real(symStream),imag(symStream),'*')
hold off
title('Received symbols VS Send symbols')
%%%

%%% Plotting the end result (Correct if it is zero everywhere)
plot_t2 = (0:numel(recStream)-1)*T;
figure()
hold on
stem(plot_t2,bitStream-recStream);
hold off
title('Original bitsream minus Recovered bitstream')
xlabel('time [s]')
ylabel('Correct if = 0')
%%%
