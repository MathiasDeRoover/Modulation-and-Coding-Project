clear
close all
clc
addpath(genpath(pwd))
%% INITIALIZATION %%
modu = 'BPSK';
[modulation,bps]    = ModuToModulation(modu);
ftaps               = 200;              % Amount of causal (and non causal) filter taps
symRate             = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T                   = 1/symRate;        % Symbol period
M                   = 100;              % UpSample factor
fs                  = symRate * M;      % Sample Frequency
beta                = 0.3;              % Roll off factor
N                   = 128*120;          % Amount of bits in original stream
hardDecodeIter      = 10;               % Iteration limit for hard decoder
cLength             = 128;
vLength             = 256;
SNRdB               = 200;              % Signal to noise in dB
fc                  = 2e9;              % 2 GHz carrier freq is given as example in the slides
ppm                 = 2e-6;             % 2 parts per million

deltaW              = 0;
phi0                = 0;                % Phase offset
delta               = 0;                % Sample clock offset SCO
t0                  = 0;                % Time shift
K_gard              = 0;                % K for Gardner


t0                  = [0 0.1 0.2 0.3 0.4]*T;

SNRdB               = -10:0.5:20;

BER                 = zeros(size(SNRdB));
%% Create bitstream
[stream_bit]        = CreateBitStream(N,1);

%% Create H matrix
H0                  = CreateHMatrix(cLength,vLength);

%% LDPC encoding
[stream_coded,newH] = LDPCencode(stream_bit,H0);

%% Mapping
stream_mapped       = mapping(stream_coded, bps, modulation);

%% Add pilot symbols

%% Upsample
stream_upSampled    = upsample(stream_mapped,M);

%% Create window
[g,g_min]           = CreateWindow(T, fs, ftaps, beta);

%% Apply window
stream_wind         = conv(stream_upSampled,g);
for tt=1:length(t0)
    for i=1:length(SNRdB)
        disp(['Current SNR: ',num2str(SNRdB(i))])
        %% Sending through channel
        noisePower          = CalcNoisePower(stream_bit, stream_wind,1/fs,SNRdB(i),fs,ftaps);
        stream_channel      = Channel(stream_wind,noisePower);

        %% CFO + Phase offset
        stream_rec_CFOPhase = AddCFOAndPhase(stream_channel,fs,deltaW,phi0);

        %% Apply window
        stream_rec_wind     = conv(stream_rec_CFOPhase,g_min);

        %% Truncate & Sample + Sample clock offset and time shift
        stream_rec_sample   = TruncateAndSample(stream_rec_wind,ftaps,T,fs,delta,t0(tt)); 



        %% Pilot ToA estimation

        %% Demapping
        stream_rec_demapped = demapping(stream_rec_sample, bps, modulation);

        %% LDPC decoding
        stream_rec_decoded  = LDPC_decoHardVec( stream_rec_demapped, newH ,hardDecodeIter);
        [~, BER(i)] = biterr(stream_rec_decoded,stream_bit);
    end
    BER_t0(tt,:) = BER;
end


figure
for i=1:length(t0)
    semilogy(SNRdB,BER_t0(i,:))
    hold on
end
title('BER in function of time shift (BPSK)')
xlabel('SNR (dB)')
ylabel('BER')
legend('t0 = 0','t0 = 0.1*T','t0 = 0.2*T','t0 = 0.3*T','t0 = 0.4*T')


rmpath(genpath(pwd))