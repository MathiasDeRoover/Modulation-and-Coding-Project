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
N                   = 128*12;            % Amount of bits in original stream
hardDecodeIter      = 10;               % Iteration limit for hard decoder
cLength             = 128;
vLength             = 256;
SNRdB               = 4;              % Signal to noise in dB
fc                  = 2e9;              % 2 GHz carrier freq is given as example in the slides
ppm                 = 2e-6;             % 2 parts per million
% deltaW              = fc*ppm;           % Carrier frequency offset CFO 10ppm 10e-6
deltaW              = 0;
phi0                = 0;                % Phase offset
delta               = 0;                % Sample clock offset SCO
t0                  = 0.5*T;                % Time shift
K_gard              = 0.01;             % K for Gardner
K_gard2              = 0.05;             % K for Gardner

% SNRdB               = -10:0.1:10;

t0 = [0 0.1 0.2 0.3 0.4]*T;

nrOfIterations = 10;                % To take mean and std
BER             = zeros(size(SNRdB));
BER_no_gardner  = zeros(size(SNRdB));
for tt=1:length(t0);
    for i=1:length(SNRdB)
        for ii=1:nrOfIterations
            %% Create bitstream
            [stream_bit]        = CreateBitStream(N,1);

            %% Create H matrix
    %         H0                  = CreateHMatrix(cLength,vLength);

            %% LDPC encoding
    %         [stream_coded,newH] = LDPCencode(stream_bit,H0);
            stream_coded = stream_bit; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            %% Mapping
            stream_mapped       = mapping(stream_coded, bps, modulation);

            %% Upsample
            stream_upSampled    = upsample(stream_mapped,M); %%!!!!!!!!!!!!!!

            %% Create window
            [g,g_min]           = CreateWindow(T, fs, ftaps, beta);

            %% Apply window
            stream_wind         = conv(stream_upSampled,g);

            %% Sending through channel
            noisePower          = CalcNoisePower(stream_bit, stream_wind,1/fs,SNRdB(i),fs,ftaps); %%%!!!!!!!!!!!!!!!
            stream_channel      = Channel(stream_wind,noisePower);
    %         stream_channel = stream_wind; %%%%!!!!!%%%%

            %% CFO + Phase offset
            stream_rec_CFOPhase = AddCFOAndPhase(stream_channel,fs,deltaW,phi0);

            %% Apply window
            stream_rec_wind     = conv(stream_rec_CFOPhase,g_min);

            %% Truncate & Sample + Sample clock offset and time shift
            stream_rec_sample   = TruncateAndSample(stream_rec_wind,ftaps,T,fs,delta,t0); %@Mathias, ik geef gwn T mee ipv T/2 om te kunnen vergelijken

            %% Gardner
    %         epsilon = zeros(nrOfIterations,numel(stream_rec_sample));        % DON'T DO THIS, THIS GIVES STRANGE RESULTS BECAUSE OF MEAN AND STD
            disp(['Gardner iteration: ',num2str(ii)])
            [stream_rec_gardner, T_est, epsilon(ii,:)]  = Gardner(stream_rec_sample, K_gard, stream_rec_wind, ftaps, fs, T,delta, t0);
    %         [~, ~, epsilon2(ii,:)]  = Gardner(stream_rec_sample, K_gard2, stream_rec_wind, ftaps, fs, T,delta, t0);
       end

    error = (t0/T-abs(epsilon))*T;
    error2 = (t0/T-abs(epsilon2))*T;
    
    avg = mean(error,1);
    stddev = std(error,1);
    avg2 = mean(error2,1);
    stddev2 = std(error2,1);
    
%     avgFilt = 1/9*ones(9,1);
%     avg     = conv(avg,avgFilt,'same');
%     stddev  = conv(stddev,avgFilt,'same');
%     avg2     = conv(avg2,avgFilt,'same');
%     stddev2  = conv(stddev2,avgFilt,'same');
    
%     avg     = avg(5:end);
%     stddev  = stddev(5:end);   
%     avg2     = avg2(5:end);
%     stddev2  = stddev2(5:end);  
    
    
    
    steps = 1:20:numel(avg);
    figure
    plot(steps,avg(steps),'-or')
    hold on
    plot(steps,avg2(steps),'-xb')
    
    plot(steps,avg(steps)-stddev(steps),'--or')
    plot(steps,avg(steps)+stddev(steps),'--or')
    plot(steps,avg2(steps)-stddev2(steps),'--xb')
    plot(steps,avg2(steps)+stddev2(steps),'--xb')
    
    title('QPSK, 2 Mbps symbol rate, 0.3 roll-off, no noise')
    xlabel('Symbols')
    ylabel('Time error (mean \pm stdv)')
    ylim(1e-7*[-0.5,3])
    legend(['K = ',num2str(K_gard)],['K = ',num2str(K_gard2)])
    hold off



    %% Demapping
    % stream_rec_demapped = demapping(stream_rec_toa, bps, modulation);
    stream_rec_demapped = demapping(stream_rec_gardner, bps, modulation);
    demapped_no_gardner = demapping(stream_rec_sample, bps, modulation);

    %% LDPC decoding
    stream_rec_decoded  = LDPC_decoHardVec( stream_rec_demapped, newH ,hardDecodeIter);
    decoded_no_gardner  = LDPC_decoHardVec( demapped_no_gardner, newH ,hardDecodeIter);
    [~, BER(i)] = biterr(stream_rec_decoded,stream_bit);
    [~, BER_no_gardner(i)] = biterr(decoded_no_gardner,stream_bit);

end




rmpath(genpath(pwd))