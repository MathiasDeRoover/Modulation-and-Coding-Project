clear
close all
clc
addpath(genpath(pwd))
%% INITIALIZATION %%
modu = 'QPSK';
[modulation,bps]    = ModuToModulation(modu);
ftaps               = 80;               % Amount of causal (and non causal) filter taps
symRate             = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T                   = 1/symRate;        % Symbol period
M                   = 10;               % UpSample factor
fs                  = symRate * M;      % Sample Frequency
beta                = 0.3;              % Roll off factor
N                   = 128*6*4;            % Amount of bits in original stream
hardDecodeIter      = 10;               % Iteration limit for hard decoder
cLength             = 128;
vLength             = 256;
% SNRdB = linspace(0,16,10);            % Signal to noise ration in dB
SNRdB = 20;

% Synchronization parameters
fc                  = 2e9;              % 2 GHz carrier freq is given as example in the slides
ppm                 = 1e-6;             % Parts per million
deltaW              = 2*fc*ppm;         % Carrier frequency offset CFO (example: 2 ppm)
phi0                = pi/2;             % Phase offset: Note that the phase offset is not compensated between the first 2 pilot sequences
delta               = 0;                % Sample clock offset SCO, should always be 0
t0                  = 0.1*T;            % Time shift
K                   = 0.05;             % K for Gardner
pilotLength         = 20; 
dataBlockLength     = 500;
K_toa               = 8;


%% Create bitstream
[stream_bit]        = CreateBitStream(N,1);

%% Create H matrix
H0                  = CreateHMatrix(cLength,vLength);

%% Create window
[g,g_min]           = CreateWindow(T, fs, ftaps, beta);

%% LDPC encoding
[stream_coded,newH] = LDPCencode(stream_bit,H0);

%% Mapping
stream_mapped       = mapping(stream_coded, bps, modulation);

for i=1:length(SNRdB)
    %% Add pilot symbols
    [stream_mapped_pilots,paddLength,pilot] = PilotAndPadd(stream_mapped,dataBlockLength,pilotLength,modu);

    %% Upsample
    stream_upSampled    = upsample(stream_mapped_pilots,M);

    %% Apply window
    stream_wind         = conv(stream_upSampled,g);

    smallResult = zeros(1,100);

    %% Sending through channel
    noisePower          = CalcNoisePower(stream_bit,stream_wind,1/fs,SNRdB(i),fs,ftaps);
    stream_channel      = Channel(stream_wind,noisePower);

    %% CFO + Phase offset
    stream_rec_CFOPhase = AddCFOAndPhase(stream_channel,fs,deltaW,phi0);

    %% Apply window
    stream_rec_wind     = conv(stream_rec_CFOPhase,g_min);

    %% Gardner
    [stream_rec_gardner, T_est]  = Gardner(stream_mapped_pilots, K, stream_rec_wind, ftaps, fs, T,delta, t0);

    %% Pilot ToA estimation
    [stream_rec_toa,deltaf,n_est] = dePilotAndDePadd(stream_rec_gardner,paddLength,pilotLength,dataBlockLength,pilot,K_toa,T_est);

    %% Demapping
    stream_rec_demapped = demapping(stream_rec_toa, bps, modulation);

    %% LDPC decoding
    stream_rec_decoded  = LDPC_decoHardVec( stream_rec_demapped, newH ,hardDecodeIter);
%     [~, BER(i)] = biterr(stream_rec_decoded,stream_bit);
end

%% Results
figure
stem(stream_bit~=stream_rec_decoded)
title('Figure of symbol errors at receiver')

% figure
% semilogy(SNRdB,BER);


rmpath(genpath(pwd))