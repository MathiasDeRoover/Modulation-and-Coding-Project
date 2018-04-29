%% Initializations
clear
close all
clc
addpath(genpath(pwd))

%% Generate H matrix
c           = 128;                    % Block length
v           = 256;                   % Coded block length
H           = makeLdpc(c,v,0,1,3);  % Create initial parity check matrix

%% Execute for multiple SNRs
SNRvec      = linspace(0,20,40);
berSoft     = zeros(size(SNRvec));
berSoftLog  = zeros(size(SNRvec));
berHard     = zeros(size(SNRvec));
ber         = zeros(size(SNRvec));

wait_bar = waitbar(0,'Please wait...');
for i = 1:numel(SNRvec)
    %% Generate bitstream
    N           = 1e3;   % Number of symbols
    bps         = 1;    % Bits per symbol
    bitStream   = CreateBitStream( N,bps );
    
    [berSoft(i),berSoftLog(i),berHard(i),ber(i)] = LDPC_berCalc(bitStream,H,SNRvec(i));
    waitbar(i/numel(SNRvec),wait_bar);
end
close(wait_bar)

figure
semilogy(SNRvec,berSoft)
hold on
semilogy(SNRvec,berSoftLog)
semilogy(SNRvec,berHard)
semilogy(SNRvec,ber)
legend('Soft Decoding','Log Soft Decoding','Hard Decoding','No Decoding');
hold off
rmpath(genpath(pwd))