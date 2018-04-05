%% Initializations
clear
close all
clc
addpath(genpath(pwd))

N        = 1e3;              % Amount of Symbols
window_N = 10;               % Time axis window = window_N * T
symRate  = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T        = 1/symRate;        % Symbol period
beta     = 0.3;              % Roll off factor
M        = [4 5 6 7 8 35 50 80];              % UpSample factors


SNRdb   = linspace(-20,20,500);
SNR     = db2mag(2*SNRdb);      %takes db:20log10(x) and we have power so it should be 10log10(x)

BPSKarray  = cell(numel(M),1);  % Initializing the data cell arrays
QPSKarray  = cell(numel(M),1);  %
Qam16array = cell(numel(M),1);  %
Qam64array = cell(numel(M),1);  %

%% Execution
start = tic;
for j=1:numel(M)
    disp(['loop: ',num2str(j), ' started'])
    tic

    fs=symRate .* M(j);
    
    %-- Filter desgin for specific M --%
    ftaps = window_N*M(j);                              % Amount of causal (and non causal) filter taps
    H = RRCFilter( T,fs,beta,ftaps )';                  % Creating the window
    h = ifft(H,'symmetric');                            % Transforming to the time domain
    h = h/h(1);                                         % Normalizing the window in the time domain
    H = fft(h);                                         % Transforming the normalized window back to the frequency domain
    
    G = sqrt(H);                                        % G is the square root of H so that G*G = H (after the convolutions)
    g = fftshift(ifft(G,'symmetric'));                  % Transform G to the time domain. Fftshift is needed to get a proper raised cosine
    
    plotFilter( fs,h,H,G,M(j),ftaps )                   % Plot H,h,G and g
    
    %-- Calculating BER --%
    BPSKarray{j}=BERcalc(N,ftaps,'BPSK',M(j),fs,g,SNR, j, numel(M));
    QPSKarray{j}=BERcalc(N,ftaps,'QPSK',M(j),fs,g,SNR, j, numel(M));
    Qam16array{j}=BERcalc(N,ftaps,'16QAM',M(j),fs,g,SNR, j, numel(M));
    Qam64array{j}=BERcalc(N,ftaps,'64QAM',M(j),fs,g,SNR, j, numel(M));
    
    toc
end
toc(start)

%%%
%Validation( shsStream,symStream,bitStream,recStream,T ) %Function to validate (plot) channel performance
%%%


%% Plotting the graphs 
figure
subplot(2,2,1)
semilogy(SNRdb,cell2mat(BPSKarray)');
xlabel('SNR=Eb/N0(dB)')
ylabel('BER')
legend(num2str(M'))
title('BER BPSK for different upsampling rates')

subplot(2,2,2)
semilogy(SNRdb,cell2mat(QPSKarray)');
xlabel('SNR=Eb/N0(dB)')
ylabel('BER')
legend(num2str(M'))
title('BER QPSK for different upsampling rates')

subplot(2,2,3)
semilogy(SNRdb,cell2mat(Qam16array)');
xlabel('SNR=Eb/N0(dB)')
ylabel('BER')
legend(num2str(M'))
title('BER Qam16 for different upsampling rates')

subplot(2,2,4)
semilogy(SNRdb,cell2mat(Qam64array)');
xlabel('SNR=Eb/N0(dB)')
ylabel('BER')
legend(num2str(M'))
title('BER Qam64 for different upsampling rates')

%suptitle('BER for different constellations and different upsampling rates')

figure
semilogy(SNRdb,BPSKarray{1},'b');
hold on
plot([4 6 8],[2e-2 3e-3 3e-4]./log2(2),'bo');
semilogy(SNRdb,QPSKarray{1},'Color',[0.80 0.40 0.0]);
plot([1 5 7],[1e-1 1e-2 1e-3]./log2(4),'o','Color',[0.80 0.40 0.0]);
semilogy(SNRdb,Qam16array{1},'r');
plot([9 10 12],[2e-2 5e-3 1e-3]./log2(16),'ro');
semilogy(SNRdb,Qam64array{1},'k');
plot([12 14 18],[5e-2 1e-2 1e-4]./log2(64),'ko');
hold off
xlabel('SNR=Eb/N0(dB)')
ylabel('BER')
legend('BPSK','BPSKth','QPSK','QPSKth','Qam16','Qam16th','Qam64','Qam64th')
title(['BER for different constellations, upsampling rate = ' num2str(M(1))])

rmpath(genpath(pwd))