%% Initializations
clear
close all
clc
addpath(genpath(pwd))

%% Generate H matrix
c_nodes     = 128;                  % Block length
v_nodes     = 256;                  % Coded block length
SNR         = linspace(-10,10,20);  % SNR in dB for IdealChannel_exec
% SNR         = db2mag(SNRdb);
H           = makeLdpc(c_nodes,v_nodes,0,1,3);  % Create initial parity check matrix

berSoft            = zeros(size(SNR));
berSoftLog         = zeros(size(SNR));
berHard            = zeros(size(SNR));
berNotCoded        = zeros(size(SNR));

%% Generate bitstream
N           = 500;   % Number of symbols
bps         = 1;    % Bits per symbol

for i=1:length(SNR)
bitStream   = CreateBitStream( N,bps );

%% Add padding if bitstream is no multiple of c
zerosToPad = c_nodes - mod(length(bitStream),c_nodes);
if (zerosToPad ~= 0)
    bitStreamPad = [bitStream; zeros(zerosToPad,1)];
    while mod(numel(bitStreamPad),bps) ~= 0                % Infinite loops should not happen, but be careful!
        bitStreamPad = [bitStreamPad; zeros(c_nodes,1)];
        zerosToPad = zerosToPad + c_nodes;
    end
else
    bitStreamPad = bitStream;
end

%% Encode bitstream
bitSBlock       = reshape(bitStreamPad,c_nodes,[]);
[bitPBlock,Hs]  = makeParityChk(bitSBlock,H,0);     % Encode & create parity check bits 
bitSCBlock      = [bitPBlock;bitSBlock];            % Concatenate parity check bits and data bits
bitSCoded       = reshape(bitSCBlock,[],1);

%% Send through channel
    receivedUncoded     = IdealChannel_exec(bitStream,SNR(i),'BPSK','no_det');
    receivedUncoded     = receivedUncoded>0;
%     receivedUncoded     = receivedUncoded(1:end-zerosToPad); No padding
%     is added for uncoded message

    receivedCoded       = IdealChannel_exec(bitSCoded,SNR(i),'BPSK','no_det');     % Send through channel withouth detector
    receivedCodedDet    = receivedCoded>0;                                           % Detector
    
%% Decode bitstream
    sigNoise            = std(receivedCoded(receivedCoded>0));
    
    tic
    bitRecoveredSoft    = LDPC_decoder_soft_BPSK( receivedCoded,Hs,sigNoise );
    bitRecoveredSoft    = bitRecoveredSoft(1:end-zerosToPad);
    toc
    
    tic
    bitRecoveredSoftLog = LDPC_decoder_soft_log_BPSK( receivedCoded,Hs,sigNoise );
    bitRecoveredSoftLog = bitRecoveredSoftLog(1:end-zerosToPad);
    toc
    
    tic
    bitRecoveredHard    = LDPC_decoder_hard( receivedCodedDet, Hs,10 );
    bitRecoveredHard    = bitRecoveredHard(1:end-zerosToPad);
    toc
    


%% Calculate bit error
    [~,berSoft(i)]     = biterr(bitRecoveredSoft,bitStream);
    [~,berSoftLog(i)]  = biterr(bitRecoveredSoftLog,bitStream);
    [~,berHard(i)]     = biterr(bitRecoveredHard,bitStream);
    [~,berNotCoded(i)] = biterr(receivedUncoded,bitStream);
    
%     berSoftM(i)     = berSoft;
%     berSoftLogM(i)  = berSoftLog;
%     berHardM(i)     = berHard;
%     berNotCodedM(i) = berNotCoded;
    
    
end

%% Plots
figure
semilogy(SNR, berNotCoded)
hold on
semilogy(SNR, berHard)
semilogy(SNR, berSoft)
semilogy(SNR, berSoftLog)
hold off
legend('Not coded','Hard','Soft','Soft Log')
xlabel('SNR (dB)')
ylabel('BER')


rmpath(genpath(pwd))