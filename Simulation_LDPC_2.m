%% Initializations
clear
close all
clc
addpath(genpath(pwd))

%% Generate H matrix
c           = 5;                  % Block length
v           = 10;                  % Coded block length
H           = makeLdpc(c,v,0,1,3);  % Create initial parity check matrix

%% Generate bitstream
N           = 50;   % Number of symbols
bps         = c;    % Bits per symbol
bitStream   = CreateBitStream( N,bps );

%% Encode bitstream
bitSBlock       = reshape(bitStream,c,[]);
[bitPBlock,Hs]  = makeParityChk(bitSBlock,H,0);     % Encode & create parity check bits 
bitSCBlock      = [bitPBlock;bitSBlock];            % Concatenate parity check bits and data bits
bitSCoded       = reshape(bitSCBlock,[],1);

%% Send through channel
SNR         = 3; 
bitSRecv    = IdealChannel_exec(bitSCoded,SNR,'BPSK','no_det');     % Send through channel withouth detector
bitSRDet    = bitSRecv>0;                                           % Detector

%% Decode bitstream
sigNoise            = std(bitSRecv(bitSRecv>0));
tic
bitRecoveredSoft    = LDPC_decoder_soft_BPSK( bitSRecv,Hs,sigNoise );
toc
tic
bitRecoveredSoftLog = LDPC_decoder_soft_log_BPSK( bitSRecv,Hs,sigNoise );
toc
tic
bitRecoveredHard    = LDPC_decoder_hard( bitSRDet, Hs,10 );
toc

%% Calculate bit error
[~,berSoft]     = biterr(bitRecoveredSoft,bitStream);
[~,berSoftLog]  = biterr(bitRecoveredSoftLog,bitStream);
[~,berHard]     = biterr(bitRecoveredHard,bitStream);

rmpath(genpath(pwd))