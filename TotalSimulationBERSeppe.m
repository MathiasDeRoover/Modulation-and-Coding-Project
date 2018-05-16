clear
close all
clc
addpath(genpath(pwd))
%% INITIALIZATION %%
modu = 'BPSK';
[modulation,bps]    = ModuToModulation(modu);
ftaps               = 800;               % Amount of causal (and non causal) filter taps
symRate             = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T                   = 1/symRate;        % Symbol period
M                   = 100;               % UpSample factor
fs                  = symRate * M;      % Sample Frequency
beta                = 0.3;              % Roll off factor
N                   = 128*6*10;              % Amount of bits in original stream
hardDecodeIter      = 10;               % Iteration limit for hard decoder
cLength             = 128;
vLength             = 256;
deltaW              = 0;                % Carrier frequency offset CFO 10ppm 10e-6
phi0                = 0;                % Phase offset
delta               = 0;                % Sample clock offset SCO
% t0                  = 0.4*T;                % Time shift
t0 = 0;% Check influence without time shift
K                   = 0.01;              % K for Gardner
fc                  = 2e9;              % 2 GHz carrier freq is given as example in the slides
pilotLength         = 20; 
dataBlockLength     = 500;
K_toa               = 10;

SNRdB               = 0:10;             % Signal to noise in dB
% SNRdB = 200;
ppm = [0 2e-6 10e-6];
% ppm = 0;
deltaW =ppm.*fc;
%% Create bitstream
[stream_bit]        = CreateBitStream(N,1);

%% Create H matrix
H0                  = CreateHMatrix(cLength,vLength);

%% LDPC encoding
[stream_coded,newH] = LDPCencode(stream_bit,H0);

%% Mapping
stream_mapped       = mapping(stream_coded, bps, modulation);

%% Add pilot symbols
[stream_mapped_pilots,paddLength,pilot] = PilotAndPadd(stream_mapped,dataBlockLength,pilotLength,modu);

%% Upsample
stream_upSampled    = upsample(stream_mapped_pilots,M);

%% Create window
[g,g_min]           = CreateWindow(T, fs, ftaps, beta);

h = conv(g,g_min,'same');
% figure
% th = 0:1/fs:(numel(h)-1)*1/fs;
% plot(th,h)
% hold on
% th2 = 0:T:(numel(h(1:T/1*fs:end))-1)*T;
% stem(th2,h(1:T/1*fs:end));
% hold off


%% Apply window
stream_wind         = conv(stream_upSampled,g);

BER = zeros(size(SNRdB));

figure

for j = 1:length(deltaW)
    disp(['CFO is: ',num2str(deltaW(j))])
    for i=1:length(SNRdB)
        %% Sending through channel
        disp(['SNR is: ',num2str(SNRdB(i))])

        noisePower          = CalcNoisePower(stream_bit,stream_wind,1/fs,SNRdB(i),fs,ftaps); %%%!!!!!!!!!!!!!!!
        stream_channel      = Channel(stream_wind,noisePower);

        %% CFO + Phase offset
        stream_rec_CFOPhase = AddCFOAndPhase(stream_channel,fs,deltaW(j),phi0);

        %% Apply window
        stream_rec_wind     = conv(stream_rec_CFOPhase,g_min);

        %% Truncate & Sample + Sample clock offset and time shift
%         stream_rec_sample   = TruncateAndSample(stream_rec_wind,ftaps,T,fs,delta,t0); % Twice as many samples are needed to implement Gardner, so we will sample here at T/2. Alternative: use upsample factor M = 2 and sample here at T.
        stream_rec_wind = stream_rec_wind(2*ftaps+1:end-2*ftaps); %No Garnder
        stream_rec_sample   = stream_rec_wind(1:M:end);
        %% Gardner
%         [stream_rec_gardner,T_est, epsilon]  = Gardner(stream_rec_sample, K, stream_rec_wind, ftaps, fs, T,delta, t0);
    %     stream_rec_gardner = stream_rec_sample(1:2:end);

        %% Plots

        plotBegin = 500;
        plotEnd   = 550;

    %     figure
    %     hold on
    %     t = 0:1/fs:(numel(stream_rec_wind(2*ftaps+1:end))-1)/fs;
    %     stream_rec_wind_plot = real(stream_rec_wind(2*ftaps+1:end));
    % %     plot(t,stream_rec_wind_plot); %full plot
    %     plot(t(plotBegin*M:plotEnd*M),stream_rec_wind_plot(plotBegin*M:plotEnd*M));
    % 
    %     t2 = 0:T/2:(numel(stream_rec_sample)-1)*T/2;
    %     stream_rec_sample_plot = real(stream_rec_sample);
    % %     stem(t2,stream_rec_sample_plot); %full plot
    %     stem(t2(plotBegin*2:plotEnd*2),stream_rec_sample_plot(plotBegin*2:plotEnd*2),'*');
    % 
    %     t3 = 0:T:(numel(stream_rec_gardner)-1)*T;
    %     stream_rec_gardner_plot = real(stream_rec_gardner);
    % %     stem(t3,stream_rec_gardner_plot); %full plot
    %     stem(t3(plotBegin:plotEnd),stream_rec_gardner_plot(plotBegin:plotEnd));
    %     legend('Convolved with matched filter','Sampled at T/2 with time shift and SCO','Reconstructed stream Gardner')
    %     hold off

    %     %Check symbol error
    %     symbolDifference = abs(stream_mapped - stream_rec_gardner);
    %     avg = 1/50*ones(1,50);
    %     symbolDifference = conv(symbolDifference,avg,'same');
    %     figure
    %     plot(symbolDifference)
    %     title('Symbol error after Gardner')

        %% Compensate only phase, keep ISI
        T_est = (0:numel(stream_mapped_pilots)-1)*T;
        [stream_rec_toa,~,~] = removePhase(stream_rec_sample,paddLength,pilotLength,dataBlockLength,pilot,K_toa,T_est,'plot');
     

        %% Demapping
    %     stream_rec_demapped = demapping(stream_rec_gardner, bps, modulation);
        demapped_no_gardner = demapping(stream_rec_toa, bps, modulation);

        %% LDPC decoding
    %     stream_rec_decoded  = LDPC_decoHardVec( stream_rec_demapped, newH ,hardDecodeIter);
        stream_rec_decoded  = LDPC_decoHardVec( demapped_no_gardner, newH ,hardDecodeIter); % No Gardner
        [~, BER(i)] = biterr(stream_rec_decoded,stream_bit);
    end

%% Results
semilogy(SNRdB,BER)
hold on
title('BER')
end
legend('CFO = 0','CFO  = 2ppm','CFO = 10ppm')
figure
stem(stream_rec_decoded~=stream_bit);
ylim([-1.1,1.1])
title('Errors after demapping and decoding')

rmpath(genpath(pwd))