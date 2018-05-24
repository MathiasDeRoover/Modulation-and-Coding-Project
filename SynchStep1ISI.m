clear
close all
clc
addpath(genpath(pwd))
%% INITIALIZATION %%
modu = '16QAM';
[modulation,bps]    = ModuToModulation(modu);
ftaps               = 200;              % Amount of causal (and non causal) filter taps
symRate             = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T                   = 1/symRate;        % Symbol period
M                   = 100;              % UpSample factor
fs                  = symRate * M;      % Sample Frequency
beta                = 0.3;              % Roll off factor
N                   = 128*500;          % Amount of bits in original stream
hardDecodeIter      = 10;               % Iteration limit for hard decoder
cLength             = 128;
vLength             = 256;
SNRdB               = 20;               % Signal to noise in dB
fc                  = 2e9;              % 2 GHz carrier freq is given as example in the slides
ppm                 = 1e-6;             % 2 parts per million
deltaW              = 2*fc*ppm;         % Carrier frequency offset CFO
phi0                = 0;                % Phase offset
delta               = 0;                % Sample clock offset SCO
t0                  = 0;                % Time shift

deltaW0 = 0;
deltaW2 = 2*fc*ppm;
deltaW10 = 10*fc*ppm;

SNRdB = -10:0.5:10;
% SNRdB = 200;


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

%% Create window
[g,g_min]           = CreateWindow(T, fs, ftaps, beta);

%% Apply window
stream_wind         = conv(stream_upSampled,g);

BER0 = zeros(size(SNRdB));
BER2 = zeros(size(SNRdB));
BER10 = zeros(size(SNRdB));

for i=1:length(SNRdB)
    %% Sending through channel
    disp(['SNRdB = ',num2str(SNRdB(i))])
    noisePower          = CalcNoisePower(stream_bit,stream_wind,1/fs,SNRdB(i),fs,ftaps);
    stream_channel      = Channel(stream_wind,noisePower);

    %% CFO + Phase offset
    stream_rec_CFOPhase0 = AddCFOAndPhase(stream_channel,fs,deltaW0,phi0);
    stream_rec_CFOPhase2 = AddCFOAndPhase(stream_channel,fs,deltaW2,phi0);
    stream_rec_CFOPhase10 = AddCFOAndPhase(stream_channel,fs,deltaW10,phi0);

    %% Apply window

    stream_rec_wind0     = conv(stream_rec_CFOPhase0,g_min);
    stream_rec_wind2     = conv(stream_rec_CFOPhase2,g_min);
    stream_rec_wind10     = conv(stream_rec_CFOPhase10,g_min);
    
    taxis = (0:numel(stream_rec_wind0)-1)*1/fs;
    stream_rec_wind0 = stream_rec_wind0 .* exp(-1i*2*pi*deltaW0*taxis.');         %% Influence of only ISI
    stream_rec_wind2 = stream_rec_wind2 .* exp(-1i*2*pi*deltaW2*taxis.');         %% Influence of only ISI
    stream_rec_wind10 = stream_rec_wind10 .* exp(-1i*2*pi*deltaW10*taxis.');      %% Influence of only ISI

    %% Truncate & Sample + Sample clock offset and time shift
    stream_rec_sample0   = TruncateAndSample(stream_rec_wind0,ftaps,T,fs,delta,t0); 
    stream_rec_sample2   = TruncateAndSample(stream_rec_wind2,ftaps,T,fs,delta,t0);
    stream_rec_sample10   = TruncateAndSample(stream_rec_wind10,ftaps,T,fs,delta,t0);
    
%     figure
%     scatter(real(stream_rec_sample2),imag(stream_rec_sample2))
%     title('16QAM symbols after adding 2 ppm CFO (noiseless)')
%     xlabel('Real part')
%     ylabel('Imaginary part')
%     defaultPos = get(0,'defaultfigureposition');
%     fig = gcf;
%     set(fig,'Position',defaultPos*2);

    %% Demapping
    stream_rec_demapped0 = demapping(stream_rec_sample0, bps, modulation);
    stream_rec_demapped2 = demapping(stream_rec_sample2, bps, modulation);
    stream_rec_demapped10 = demapping(stream_rec_sample10, bps, modulation);

    %% LDPC decoding
    stream_rec_decoded0  = LDPC_decoHardVec( stream_rec_demapped0, newH ,hardDecodeIter);
    stream_rec_decoded2  = LDPC_decoHardVec( stream_rec_demapped2, newH ,hardDecodeIter);
    stream_rec_decoded10  = LDPC_decoHardVec( stream_rec_demapped10, newH ,hardDecodeIter);
    
    [~, BER0(i)] = biterr(stream_rec_decoded0,stream_bit);
    [~, BER2(i)] = biterr(stream_rec_decoded2,stream_bit);
    [~, BER10(i)] = biterr(stream_rec_decoded10,stream_bit);
end
%% Results

figure
semilogy(SNRdB,BER0,'b')
hold on
semilogy(SNRdB,BER2,'rx')
semilogy(SNRdB,BER10,'ko')
hold off
legend('no cfo','2 ppm cfo','10 ppm cfo')
title('Impact of CFO on BER for 16QAM (only ISI)')
xlabel('SNR (dB)')
ylabel('BER')
grid on

% figure
% plot(stream_rec_decoded ~= stream_bit)

rmpath(genpath(pwd))