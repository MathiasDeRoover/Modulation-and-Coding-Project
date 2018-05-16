clear
close all
clc
addpath(genpath(pwd))
%% INITIALIZATION %%
modu = 'QPSK';
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
SNRdB               = 200;              % Signal to noise in dB
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
nrOfIterations = 10;                % To take mean and std
BER             = zeros(size(SNRdB));
BER_no_gardner  = zeros(size(SNRdB));
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

        %% Add pilot symbols
        % pilotL = 24;            % Note: this value should be divisable by 6 in case of 64QAM
        % % pilotFreq = 0.10;
        % nrOfPilots = 3;
        % [stream_mapped_pilots, pilot] = addPilots(stream_mapped, pilotL, nrOfPilots, 'BPSK');

        % Test
        % testPilots = zeros(size(stream_mapped));
        % [testPlot, ~] = addPilots(testPilots,pilotL, nrOfPilots, modu);
        % figure
        % stem(abs(testPlot))
        % ylim([-1.1,1.1])

        %% Upsample
        stream_upSampled    = upsample(stream_mapped,M); %%!!!!!!!!!!!!!!

        %% Create window
        [g,g_min]           = CreateWindow(T, fs, ftaps, beta);

        %% Apply window
        stream_wind         = conv(stream_upSampled,g);

        %% Sending through channel
%         noisePower          = CalcNoisePower(stream_bit, stream_wind,1/fs,SNRdB(i),fs,ftaps); %%%!!!!!!!!!!!!!!!
%         stream_channel      = Channel(stream_wind,noisePower);
        stream_channel = stream_wind; %%%%!!!!!%%%%

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
        [~, ~, epsilon2(ii,:)]  = Gardner(stream_rec_sample, K_gard2, stream_rec_wind, ftaps, fs, T,delta, t0);
   end

    error = (t0/T-abs(epsilon))*T;
    error2 = (t0/T-abs(epsilon2))*T;
    
    avg = mean(error,1);
    stddev = std(error,1);
    avg2 = mean(error2,1);
    stddev2 = std(error2,1);
    
    avgFilt = 1/9*ones(9,1);
    avg     = conv(avg,avgFilt,'same');
    stddev  = conv(stddev,avgFilt,'same');
    avg2     = conv(avg2,avgFilt,'same');
    stddev2  = conv(stddev2,avgFilt,'same');
    
    avg     = avg(5:end);
    stddev  = stddev(5:end);   
    avg2     = avg2(5:end);
    stddev2  = stddev2(5:end);  
    
    
    
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

    pause

    %% Pilot ToA estimation
    % K = 10;
    % k = 1:K;
    % k = 1;
    % I_n = stream_rec_gardner;
    % Dk = zeros(size(I_n));
    % % 
    % % % N in slides is the lenght of the pilot sequence
    % % 
    % Cte = 1/(pilotL-k)*exp(-1j*deltaW*k*T);
    % tic
    % for n=1:length(I_n)-pilotL                  % -pilotL such that I_n(n+l+1) cannot go out of bounds
    %     Dk_sum = zeros(length(k:pilotL-1),1);
    %     for l=k:pilotL-1
    %         Dk_sum(l) = conj(I_n(n+l+1))*I_n(n+l-k+1)*pilot(l+1)*conj(pilot(l-k+1)); % +1 indices can't be 0 in matlab
    %     end
    %     Dk(n) = Cte*sum(Dk_sum);
    % end
    % toc
    % figure
    % plot(real(Dk))
    % ylim([-1.1,1.1])
    % 
    % % Find peak indices
    % peakIndices = find(Dk>0.95);

    % Some more attempts to find peaks
    % [pks,locs] = findpeaks(real(Dk(peakIndices)))      % This never finds first peak?
    % peakIndices2 = peakIndices(locs)
    % Dk0 = zeros(size(Dk));
    % Dk0(peakIndices) = Dk(peakIndices);
    % figure
    % plot(real(Dk0))
    % title('Nice peak plot')
    % [pks,locs] = findpeaks(real(Dk0))

    %% Demapping
    % stream_rec_demapped = demapping(stream_rec_toa, bps, modulation);
    stream_rec_demapped = demapping(stream_rec_gardner, bps, modulation);
    demapped_no_gardner = demapping(stream_rec_sample, bps, modulation);

    %% LDPC decoding
    stream_rec_decoded  = LDPC_decoHardVec( stream_rec_demapped, newH ,hardDecodeIter);
    decoded_no_gardner  = LDPC_decoHardVec( demapped_no_gardner, newH ,hardDecodeIter);
    [~, BER(i)] = biterr(stream_rec_decoded,stream_bit);
    [~, BER_no_gardner(i)] = biterr(decoded_no_gardner,stream_bit);
    
%     waitbar(i/numel(SNRdB),wait_bar);
end
% close(wait_bar)

figure
semilogy(SNRdB,BER)
hold on
semilogy(SNRdB,BER_no_gardner)
hold off
legend('Gardner','No Gardner')




%% Results
% Plot filter
% h = conv(g,g_min,'same');
% figure
% th = 0:1/fs:(numel(h)-1)*1/fs;
% plot(th,h)
% hold on
% th2 = 0:T:(numel(h(1:T/1*fs:end))-1)*T;
% stem(th2,h(1:T/1*fs:end));
% hold off
% 
% Plot result
% figure
% stem(stream_rec_decoded~=stream_bit);
% ylim([-1.1,1.1])
% title('Errors after demapping and decoding')
% 
% % If you get 'Index exceeds matrix dimensions', you probably have to change
% % plotBegin and -End
% plotBegin = 500;
% plotEnd   = 520;
% 
% figure
% hold on
% time1 = 0:1/fs:(numel(stream_rec_wind(2*ftaps:end-2*ftaps+1))-1)/fs;
% stream_rec_wind_plot = real(stream_rec_wind(2*ftaps:end-2*ftaps-1));
% % plot(time1,stream_rec_wind_plot); %full plot
% plot(time1(plotBegin*M+1:plotEnd*M),stream_rec_wind_plot(plotBegin*M+1:plotEnd*M));
% 
% time2 = (0:T:(numel(stream_rec_sample)-1)*T)+t0;
% stream_rec_sample_plot = real(stream_rec_sample);
% % stem(time2,stream_rec_sample_plot); %full plot
% stem(time2(plotBegin+1:plotEnd),stream_rec_sample_plot(plotBegin+1:plotEnd),'*');
% 
% time3 = 0:T:(numel(stream_rec_gardner)-1)*T;
% stream_rec_gardner_plot = real(stream_rec_gardner);
% % stem(time3,stream_rec_gardner_plot); %full plot
% stem(time3(plotBegin+1:plotEnd),stream_rec_gardner_plot(plotBegin+1:plotEnd));
% legend('Convolved with matched filter','Sampled at T/2 with time shift','Reconstructed stream Gardner')
% 
% % plot lines at -1 and 1
% plot(time3(plotBegin+1:plotEnd),ones(size(time3(plotBegin+1:plotEnd))),'--')
% plot(time3(plotBegin+1:plotEnd),-ones(size(time3(plotBegin+1:plotEnd))),'--')
% hold off


% % Check symbol error
% symbolDifference = abs(stream_mapped - stream_rec_gardner);
% avg = 1/50*ones(1,50);
% symbolDifference = conv(symbolDifference,avg,'same');
% figure
% plot(symbolDifference)
% title('Symbol error after Gardner')






% avg = 1/50*ones(1,50);
% errorPlot = conv(e,avg,'same');
% errorPlot = e*periodToSample-timeShift;
% figure
% plot(errorPlot)
% title('Error in function of iterations')
% To remove monotonic decrease of error in case CFO and timeShift are
% perfectly known
% diff_e = diff(e);
% figure
% plot(diff_e)
% title('diff_e')

% inStreamSampled = inStream(1:2:end);
% figure
% plot(abs(inStreamSampled-outStream))

rmpath(genpath(pwd))