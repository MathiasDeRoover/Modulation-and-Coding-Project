function [berSoft,berSoftLog,berHard,ber] = LDPC_berCalc(bitStream,H,SNR)
%% Encode bitstream
bitSBlock       = reshape(bitStream,size(H,1),[]);
[bitPBlock,Hs]  = makeParityChk(bitSBlock,H,0);     % Encode & create parity check bits 
bitSCBlock      = [bitPBlock;bitSBlock];            % Concatenate parity check bits and data bits
bitSCoded       = reshape(bitSCBlock,[],1);

%% Send through channel
bitSRecv    = IdealChannel_exec(bitSCoded,SNR,'BPSK','no_det');     % Send through channel without detector
bitSRDet    = bitSRecv>0;                                           % Detector

%% Decode bitstream
sigNoise            = std(bitSRecv(bitSRecv>0));
bitRecoveredSoft    = LDPC_decoder_soft_BPSK( real(bitSRecv),Hs,sigNoise );
bitRecoveredSoftLog = LDPC_decoder_soft_log_BPSK( real(bitSRecv),Hs,sigNoise );
bitRecoveredHard    = LDPC_decoder_hard( bitSRDet, Hs,10 );

%% Calculate bit error
[~,berSoft]     = biterr(bitRecoveredSoft,bitStream);
[~,berSoftLog]  = biterr(bitRecoveredSoftLog,bitStream);
[~,berHard]     = biterr(bitRecoveredHard,bitStream);
[~,ber]         = biterr(bitSRDet,bitSCoded);
end

