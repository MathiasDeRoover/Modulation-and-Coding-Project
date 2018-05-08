clear
close all
clc
addpath(genpath(pwd))
%% INITIALIZATION %%
modu = 'BPSK';
[modulation,bps]    = ModuToModulation(modu);
ftaps               = 80;               % Amount of causal (and non causal) filter taps
symRate             = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T                   = 1/symRate;        % Symbol period
M                   = 10;               % UpSample factor
fs                  = symRate * M;      % Sample Frequency
beta                = 0.3;              % Roll off factor
N                   = 1280;              % Amount of bits in original stream
hardDecodeIter      = 10;               % Iteration limit for hard decoder
cLength             = 128;
vLength             = 256;
SNRdB               = 2000;             % Signal to noise in dB
deltaW              = 0.3*T;                % Carrier frequency offset CFO 10ppm 10e-6
phi0                = 0;                % Phase offset
delta               = 0.001;                % Sample clock offset SCO
t0                  = 0;                % Time shift
K                   = 0.1;              % K for Gardner

%% Create bitstream
[stream_bit]        = CreateBitStream(N,1);

%% Create H matrix
H0                  = CreateHMatrix(cLength,vLength);

%% LDPC encoding
[stream_coded,newH] = LDPCencode(stream_bit,H0);

%% Mapping
stream_mapped       = mapping(stream_coded, bps, modulation);

%% Upsample
stream_upSampled    = upsample(stream_mapped,M);

% IGNORE THIS
% 
% time = ((1:numel(stream_upSampled))-1)/fs;
% time2 = ((1:numel(stream_mapped))-1)*T;
% stream_downSampled = interp1(time,stream_upSampled,time2).';
% figure
% stem(time(1:40),real(stream_upSampled(1:40)))
% ylim([-1.1,1.1])
%% Create window
[g,g_min]           = CreateWindow(T/2, fs, ftaps, beta);

h = conv(g,g_min,'same');
figure
th = 0:1/fs:(numel(h)-1)*1/fs;
plot(th,h)
hold on
th2 = 0:T:(numel(h(1:T/1*fs:end))-1)*T;
stem(th2,h(1:T/1*fs:end));
hold off


%% Apply window
stream_wind         = conv(stream_upSampled,g);

%% Sending through channel
noisePower          = CalcNoisePower(stream_bit,1/fs,SNRdB,fs); %%%!!!!!!!!!!!!!!!
stream_channel      = Channel(stream_wind,noisePower);

%% CFO + Phase offset
stream_rec_CFOPhase = AddCFOAndPhase(stream_channel,fs,deltaW,phi0);
%  sum(stream_channel~=stream_rec_CFOPhase)        % if equal to 0 no CFO was added
%% Apply window
stream_rec_wind     = conv(stream_rec_CFOPhase,g_min);

%% Truncate & Sample + Sample clock offset and time shift
stream_rec_sample   = TruncateAndSample(stream_rec_wind,ftaps,T/2,fs,delta,t0); % Twice as many samples are needed to implement Gardner, so we will sample here at T/2. Alternative: use upsample factor M = 2 and sample here at T.
% stream_rec_sample   = TruncateAndSample(stream_rec_wind,ftaps,T,fs,delta,t0);

%% Gardner
stream_rec_gardner  = Gardner(stream_rec_sample, K, stream_rec_wind, ftaps, fs, T,delta, t0);

plotLim = 10;

figure
hold on
t = 0:1/fs:(numel(stream_rec_wind(2*ftaps+1:end))-1)/fs;
stream_rec_wind_plot = real(stream_rec_wind(2*ftaps+1:end));
% plot(t,stream_rec_wind_plot); %full plot
stem(t(1:plotLim*M),stream_rec_wind_plot(1:plotLim*M));

t2 = 0:T/2:(numel(stream_rec_sample)-1)*T/2;
stream_rec_sample_plot = real(stream_rec_sample);
% stem(t2,stream_rec_sample_plot); %full plot
stem(t2(1:plotLim*2),stream_rec_sample_plot(1:plotLim*2),'*');

t3 = 0:T:(numel(stream_rec_gardner)-1)*T;
stream_rec_gardner_plot = real(stream_rec_gardner);
% stem(t3,stream_rec_gardner_plot); %full plot
stem(t3(1:plotLim),stream_rec_gardner_plot(1:plotLim));
legend('Convolved with matched filter','Sampled at T/2 with time shift and SCO','Reconstructed stream Gardner')
hold off

% Check symbol error
symbolDifference = abs(stream_mapped - stream_rec_gardner);
avg = 1/50*ones(1,50);
symbolDifference = conv(symbolDifference,avg,'same');
figure
plot(symbolDifference)
title('Symbol error after Gardner (and after averaging)')
% Boundaries cause some errors
%% Demapping
stream_rec_demapped = demapping(stream_rec_gardner, bps, modulation);

%% LDPC decoding
stream_rec_decoded  = LDPC_decoHardVec( stream_rec_demapped, newH ,hardDecodeIter);

%% Results
figure
stem(stream_rec_decoded~=stream_bit);
ylim([-1.1,1.1])
title('Errors after demapping and decoding')

rmpath(genpath(pwd))