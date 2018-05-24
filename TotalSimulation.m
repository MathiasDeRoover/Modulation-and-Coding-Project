clear
% close all
clc
addpath(genpath(pwd))
%% INITIALIZATION %%
modu = 'QPSK';
[modulation,bps]    = ModuToModulation(modu);
ftaps               = 200;               % Amount of causal (and non causal) filter taps
symRate             = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T                   = 1/symRate;        % Symbol period
M                   = 10;               % UpSample factor
fs                  = symRate * M;      % Sample Frequency
beta                = 0.3;              % Roll off factor
N                   = 128*6*10;          % Amount of bits in original stream
hardDecodeIter      = 10;               % Iteration limit for hard decoder
cLength             = 128;
vLength             = 256;
SNRdB               = 20;               % Signal to noise in dB
fc                  = 2e9;              % 2 GHz carrier freq is given as example in the slides
ppm                 = 1e-6;             % 2 parts per million
deltaW              = 400*fc*ppm;           % Carrier frequency offset CFO 10ppm 10e-6  fc*ppm
phi0                = 0;                % Phase offset
delta               = 0;                % Sample clock offset SCO 0.01
t0                  = 0;                % Time shift
K                   = 0.06;             % K for Gardner
pilotLength         = 40; 
dataBlockLength     = 500;
K_toa               = 10;

%% Create bitstream
[stream_bit]        = CreateBitStream(N,1);

%% Create H matrix
H0                  = CreateHMatrix(cLength,vLength);

%% LDPC encoding
[stream_coded,newH] = LDPCencode(stream_bit,H0);

%% Mapping
stream_mapped       = mapping(stream_coded, bps, modulation);

%% Add pilot symbols
%[stream_mapped_pilots,paddLength,pilot] = PilotAndPadd(stream_mapped,dataBlockLength,pilotLength,modu);
stream_mapped_pilots=stream_mapped;
%% Upsample
stream_upSampled    = upsample(stream_mapped_pilots,M);

%% Create window
[g,g_min]           = CreateWindow(T, fs, ftaps, beta);

%% Apply window
stream_wind         = conv(stream_upSampled,g);

%% Sending through channel
noisePower          = CalcNoisePower(stream_bit,stream_wind,1/fs,SNRdB,fs,ftaps); %%%!!!!!!!!!!!!!!!
stream_channel      = Channel(stream_wind,noisePower);

%% CFO + Phase offset
stream_rec_CFOPhase = AddCFOAndPhase(stream_channel,fs,deltaW,phi0);

%% Apply window

stream_rec_wind     = conv(stream_rec_CFOPhase,g_min);
taxis = (0:numel(stream_rec_wind)-1)*1/fs;
stream_rec_wind = stream_rec_wind .* exp(-1i*2*pi*deltaW*taxis.');      %% Influence of only ISI

% For only phase drift: add cfo after 2nd window
%% Truncate & Sample + Sample clock offset and time shift
stream_rec_sample   = TruncateAndSample(stream_rec_wind,ftaps,T,fs,delta,t0); % Twice as many samples are needed to implement Gardner, so we will sample here at T/2. Alternative: use upsample factor M = 2 and sample here at T.
figure
scatter(real(stream_rec_sample),imag(stream_rec_sample))
%% Gardner
[stream_rec_gardner, T_est]  = Gardner(stream_rec_sample, K, stream_rec_wind, ftaps, fs, T,delta, t0);

%% Pilot ToA estimation
%[stream_rec_toa,~,~] = dePilotAndDePadd(stream_rec_gardner,paddLength,pilotLength,dataBlockLength,pilot,K_toa,T_est,'plot');
stream_rec_toa = stream_rec_gardner;
%% Demapping
stream_rec_demapped = demapping(stream_rec_toa, bps, modulation);

%% LDPC decoding
stream_rec_decoded  = LDPC_decoHardVec( stream_rec_demapped, newH ,hardDecodeIter);

%% Results

hold on
tmp = zeros(size(stream_mapped_pilots));
tmp(dataBlockLength+1:dataBlockLength+pilotLength:numel(stream_mapped_pilots)) = 1;
stem(tmp)
stem(linspace(1,numel(tmp),numel(stream_coded)),-(stream_rec_demapped~=stream_coded));
hold off
%ylim([-1.1,1.1])
%title('Errors after demapping and decoding')

figure
hold on
scatter(real(stream_rec_sample),imag(stream_rec_sample));
scatter(real(stream_rec_gardner),imag(stream_rec_gardner));
title('Gardner')
legend('wrong sampling','after Gardner')
hold off

figure
hold on
scatter(real(stream_rec_gardner),imag(stream_rec_gardner));
scatter(real(stream_rec_toa),imag(stream_rec_toa));
title('Pilots');
legend('after Gardner','after ToA')
hold off

rmpath(genpath(pwd))