%% Initializations
clear
close all
clc
addpath(genpath(pwd))

N        = 1024;              % Amount of Symbols
window_N = 10;               % Time axis window = window_N * T
symRate  = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T        = 1/symRate;        % Symbol period
beta     = 0.3;              % Roll off factor
M = 4;

% Limit for the nr of iterations:
% iterations=[1 5 20 100 2000 20000];
iterations=1;

SNRdb   = linspace(-20,20,500);
% SNRdb   = 60;
SNR     = db2mag(2*SNRdb);      %takes db:20log10(x) and we have power so it should be 10log10(x)

BPSKarray  = cell(numel(M),1);  % Initializing the data cell arrays
QPSKarray  = cell(numel(M),1);  %
Qam16array = cell(numel(M),1);  %
Qam64array = cell(numel(M),1);  %

%% Execution
j=1;
disp(['loop: ',num2str(j), ' started'])
tic

fs=symRate .* M;

%-- Filter design for specific M --%
ftaps = window_N*M(j);                              % Amount of causal (and non causal) filter taps
H = RRCFilter( T,fs,beta,ftaps )';                  % Creating the window
h = ifft(H,'symmetric');                            % Transforming to the time domain
h = h/h(1);                                         % Normalizing the window in the time domain
H = fft(h);                                         % Transforming the normalized window back to the frequency domain

G = sqrt(H);                                        % G is the square root of H so that G*G = H (after the convolutions)
g = fftshift(ifft(G,'symmetric'));                  % Transform G to the time domain. Fftshift is needed to get a proper raised cosine

% plotFilter( fs,h,H,G,M,ftaps )                   % Plot H,h,G and g

%-- Calculating BER without LDPC --%
BPSKarray=BERcalc(N,ftaps,'BPSK',M,fs,g,SNR, 1, numel(M));
QPSKarray=BERcalc(N,ftaps,'QPSK',M,fs,g,SNR, 1, numel(M));
Qam16array=BERcalc(N,ftaps,'16QAM',M,fs,g,SNR, 1, numel(M));
Qam64array=BERcalc(N,ftaps,'64QAM',M,fs,g,SNR, 1, numel(M));

for j=1:numel(iterations)
    BPSKarray_LDPC(j,:)=BERcalc_LDPC(N,ftaps,'BPSK',M,fs,g,SNR,j);
    QPSKarray_LDPC(j,:)=BERcalc_LDPC(N,ftaps,'QPSK',M,fs,g,SNR,j);
    Qam16array_LDPC(j,:)=BERcalc_LDPC(N,ftaps,'16QAM',M,fs,g,SNR,j);
    Qam64array_LDPC(j,:)=BERcalc_LDPC(N,ftaps,'64QAM',M,fs,g,SNR,j);
end
toc


%%%
%Validation( shsStream,symStream,bitStream,recStream,T ) %Function to validate (plot) channel performance
%%%


%% Plotting the graphs 

figure
semilogy(SNRdb,BPSKarray,'b');
hold on
semilogy(SNRdb,QPSKarray,'Color',[0.80 0.40 0.0]);
semilogy(SNRdb,Qam16array,'r');
semilogy(SNRdb,Qam64array,'k');
hold off
xlabel('SNR=Eb/N0(dB)')
ylabel('BER')
legend('BPSK','QPSK','Qam16','Qam64')
title(['BER for different constallations and M= ' num2str(M) ' without LDPC'])

figure
semilogy(SNRdb,BPSKarray_LDPC(1,:),'b');
hold on
semilogy(SNRdb,QPSKarray_LDPC(1,:),'Color',[0.80 0.40 0.0]);
semilogy(SNRdb,Qam16array_LDPC(1,:),'r');
semilogy(SNRdb,Qam64array_LDPC(1,:),'k');
hold off
xlabel('SNR=Eb/N0(dB)')
ylabel('BER')
legend('BPSK','QPSK','Qam16','Qam64')
title(['BER for different constallations and M = ' num2str(M) ' with LDPC, hard decoding'])


figure
subplot(2,2,1)
semilogy(SNRdb,BPSKarray);
hold on
semilogy(SNRdb,BPSKarray_LDPC(1,:));
legend('No LDPC','Hard decoding')
title('BPSK')
subplot(2,2,2)
semilogy(SNRdb,QPSKarray);
hold on
semilogy(SNRdb,QPSKarray_LDPC(1,:));
legend('No LDPC','Hard decoding')
title('QPSK')
subplot(2,2,3)
semilogy(SNRdb,Qam16array);
hold on
semilogy(SNRdb,Qam16array_LDPC(1,:));
legend('No LDPC','Hard decoding')
title('16QAM')
subplot(2,2,4)
semilogy(SNRdb,Qam64array);
hold on
semilogy(SNRdb,Qam64array_LDPC(1,:));
legend('No LDPC','Hard decoding')
title('64QAM')

figure
subplot(2,2,1)
legent=[]; %Legend entries
for i=1:numel(iterations)
    semilogy(SNRdb,BPSKarray_LDPC(i,:));
    hold on
    legent{i}=[num2str(iterations(i))];
end  
hold off
title('BPSK for increasing iteration limit')
legend(legent)
xlabel('SNR=E_b/N_0(dB)')
ylabel('BER')

subplot(2,2,2)
legent=[] ;
for i=1:numel(iterations)
    semilogy(SNRdb,QPSKarray_LDPC(i,:));
    hold on
    legent{i}=[num2str(iterations(i))];
end
legend(legent)
title('QPSK for increasing iteration limit')
hold off
xlabel('SNR=E_b/N_0(dB)')
ylabel('BER')

subplot(2,2,3)
legent=[] ;
for i=1:numel(iterations)
    semilogy(SNRdb,Qam16array_LDPC(i,:));
    hold on
    legent{i}=[num2str(iterations(i))];
end
legend(legent)
title('16QAM for increasing iteration limit')
hold off
xlabel('SNR=E_b/N_0(dB)')
ylabel('BER')

subplot(2,2,4)
legent=[] ;
for i=1:numel(iterations)
    semilogy(SNRdb,Qam64array_LDPC(i,:));
    hold on
    legent{i}=[num2str(iterations(i))];
end
legend(legent)
title('64QAM for increasing iteration limit')
hold off
xlabel('SNR=E_b/N_0(dB)')
ylabel('BER')

rmpath(genpath(pwd))