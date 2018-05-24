clear
close all
clc
addpath(genpath(pwd))

%% Syncronization Step4:
% The fourth step is to write a new function implementing the frame and frequency acquisition.
% Illustrate the remaining time/frequency error variances at the desired SNR for a
% varying pilot sequence length and parameter K. The robustness of the frame acquisi-
% tion to the CFO should also be checked.

modu = 'QPSK';
[modulation,bps]    = ModuToModulation(modu);
ftaps               = 80;               % Amount of causal (and non causal) filter taps
symRate             = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T                   = 1/symRate;        % Symbol period
M                   = 10;               % UpSample factor
fs                  = symRate * M;      % Sample Frequency
beta                = 0.3;              % Roll off factor
N                   = 128*6*4;          % Amount of bits in original stream
hardDecodeIter      = 10;               % Iteration limit for hard decoder
cLength             = 128;
vLength             = 256;
SNRdB               = 20;               % Signal to noise in dB
fc                  = 2e9;              % 2 GHz carrier freq is given as example in the slides
ppm                 = 1e-6;             % 2 parts per million
deltaW              = 0*fc*ppm;         % Carrier frequency offset CFO 10ppm 10e-6  fc*ppm
phi0                = 0;                % Phase offset
delta               = 0;                % Sample clock offset SCO 0.01
t0                  = 0;                % Time shift
K                   = 0.06;             % K for Gardner
pilotLength         = 20; 
dataBlockLength     = 500;
K_toa               = 8;




%% Create H matrix
H0                  = CreateHMatrix(cLength,vLength);

%% Create window
[g,g_min]           = CreateWindow(T, fs, ftaps, beta);

figure()
changeVar_vec = [0,10];
SNR_vec = linspace(0,16,10);
result = zeros(numel(changeVar_vec),numel(SNR_vec));
for j = 1:numel(changeVar_vec)
    deltaW = changeVar_vec(j)*fc*ppm;
    for i = 1:numel(SNR_vec)
        SNRdB = SNR_vec(i);
        %% Create bitstream
        [stream_bit]        = CreateBitStream(N,1);
        
        %% LDPC encoding
        [stream_coded,newH] = LDPCencode(stream_bit,H0);
        
        %% Mapping
        stream_mapped       = mapping(stream_coded, bps, modulation);
        
        %% Add pilot symbols
        [stream_mapped_pilots,paddLength,pilot] = PilotAndPadd(stream_mapped,dataBlockLength,pilotLength,modu);
        
        %% Upsample
        stream_upSampled    = upsample(stream_mapped_pilots,M);
        
        %% Apply window
        stream_wind         = conv(stream_upSampled,g);
        
        smallResult = zeros(1,100);
        for q = 1:40
            
            %% Sending through channel
            noisePower          = CalcNoisePower(stream_bit,stream_wind,1/fs,SNRdB,fs,ftaps); %%%!!!!!!!!!!!!!!!
            stream_channel      = Channel(stream_wind,noisePower);
            
            %% CFO + Phase offset
            stream_rec_CFOPhase = AddCFOAndPhase(stream_channel,fs,deltaW,phi0);
            
            %% Apply window
            stream_rec_wind     = conv(stream_rec_CFOPhase,g_min);
            
            %% Gardner
            %[stream_rec_gardner, T_est]  = Gardner(stream_mapped_pilots, K, stream_rec_wind, ftaps, fs, T,delta, t0);
            tmp = stream_rec_wind(2*ftaps-1:end-2*ftaps);
            tmp = tmp(1:M:end);
            T_est = (0:numel(stream_mapped_pilots)-1)*T;
            stream_rec_gardner = tmp;
            
            %% Pilot ToA estimation
            n_actual = dataBlockLength+1:dataBlockLength+pilotLength:numel(stream_mapped_pilots);
            [stream_rec_toa,deltaf,n_est] = dePilotAndDePadd(stream_rec_gardner,paddLength,pilotLength,dataBlockLength,pilot,K_toa,T_est);
            
            %smallResult(q) = (deltaW-deltaf)/(fc*ppm);
            smallResult(q) = std(n_actual-n_est);
        end
        
        result(j,i) = mean(smallResult);
        
        % %% Demapping
        % stream_rec_demapped = demapping(stream_rec_toa, bps, modulation);
        %
        % %% LDPC decoding
        % stream_rec_decoded  = LDPC_decoHardVec( stream_rec_demapped, newH ,hardDecodeIter);
        disp(i)
    end
end
plot(SNR_vec,result);   
xlabel('Eb/N0 [dB]')
%ylabel('Frequency error stdev [ppm]')
ylabel('Time error stdev [sampe]')
title('Random QPSK pilot symbols, no CFO, no time shift')
xlim([0,16]);


rmpath(genpath(pwd))
