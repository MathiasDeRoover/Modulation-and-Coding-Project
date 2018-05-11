%% Initializations
clear
close all
clc
addpath(genpath(pwd))

%% Generate H matrix
c_nodes     = 128;                  % Block length
v_nodes     = 256;                  % Coded block length
SNR         = linspace(-20,20,40);  % SNR in dB for IdealChannel_exec
% SNR         = db2mag(SNRdb);
H           = makeLdpc(c_nodes,v_nodes,0,1,3);  % Create initial parity check matrix

BPSK_JMG        = zeros(size(SNR));
BPSK_Guylian    = zeros(size(SNR));
UncodedBER      = zeros(size(SNR));


%% Generate bitstream
N           = 768*20;   % Number of symbols
bps         = 1;    % Bits per symbol

wait_bar = waitbar(0,'Simulating');
itlim = [1 2 5 7 9 10];
% itlim = [1 10 50 100];
for j = 1:numel(itlim)
    for i=1:length(SNR)
        bitStream   = CreateBitStream( N,bps );

    %% Add padding if bitstream is no multiple of c
    % zerosToPad = c_nodes - mod(length(bitStream),c_nodes);
    % if (zerosToPad ~= 0)
    %     bitStreamPad = [bitStream; zeros(zerosToPad,1)];
    %     while mod(numel(bitStreamPad),bps) ~= 0                % Infinite loops should not happen, but be careful!
    %         bitStreamPad = [bitStreamPad; zeros(c_nodes,1)];
    %         zerosToPad = zerosToPad + c_nodes;
    %     end
    % else
    %     bitStreamPad = bitStream;
    % end
        bitStreamPad = bitStream;

    %% Encode bitstream
        bitSBlock       = reshape(bitStreamPad,c_nodes,[]);
        [bitPBlock,Hs]  = makeParityChk(bitSBlock,H,0);     % Encode & create parity check bits 
        bitSCBlock      = [bitPBlock;bitSBlock];            % Concatenate parity check bits and data bits
        bitSCoded       = reshape(bitSCBlock,[],1);

    %% Send over channel
        receivedStream1 = IdealChannel_exec(bitSCoded,SNR(i),'BPSK','det');
        receivedStream2 = IdealChannel_exec(bitSCoded,SNR(i),'QPSK','det');
        receivedStream3 = IdealChannel_exec(bitSCoded,SNR(i),'16QAM','det');
        receivedStream4 = IdealChannel_exec(bitSCoded,SNR(i),'64QAM','det');

    %% Decode bitstream
        hardDecoded1    = LDPC_decoder_hard_biased( receivedStream1, Hs, itlim(j) );
        hardDecoded2    = LDPC_decoder_hard_biased( receivedStream2, Hs, itlim(j) );
        hardDecoded3    = LDPC_decoder_hard_biased( receivedStream3, Hs, itlim(j) );
        hardDecoded4    = LDPC_decoder_hard_biased( receivedStream4, Hs, itlim(j) );

    %% Calculate bit error
        [~,BPSK_JMG(j,i)]     = biterr(hardDecoded1,bitStream);
        [~,QPSK_JMG(j,i)]     = biterr(hardDecoded2,bitStream);
        [~,QAM16_JMG(j,i)]     = biterr(hardDecoded3,bitStream);
        [~,QAM64_JMG(j,i)]     = biterr(hardDecoded4,bitStream);
        
        waitbar(((j-1)*numel(SNR)+i)/numel(SNR)/numel(itlim),wait_bar);
    end
end
close(wait_bar);

%% Plots

% BPSK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
legent=[] ;
for k=1:numel(itlim)
semilogy(SNR, BPSK_JMG(k,:));
hold on
legent{k}=[num2str(itlim(k))];
end
hold off
legend(legent)
xlabel('SNR (dB)')
ylabel('BER')
title({'Comparison of BER for LDPC for different iteration limits','BPSK'})

% QPSK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
legent=[] ;
for k=1:numel(itlim)
semilogy(SNR, QPSK_JMG(k,:));
hold on
legent{k}=[num2str(itlim(k))];
end
hold off
legend(legent)
xlabel('SNR (dB)')
ylabel('BER')
title({'Comparison of BER for LDPC for different iteration limits','QPSK'})

% 16QAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
legent=[] ;
for k=1:numel(itlim)
semilogy(SNR, QAM16_JMG(k,:));
hold on
legent{k}=[num2str(itlim(k))];
end
hold off
legend(legent)
xlabel('SNR (dB)')
ylabel('BER')
title({'Comparison of BER for LDPC for different iteration limits','QAM16'})

% 64QAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
legent=[] ;
for k=1:numel(itlim)
semilogy(SNR, QAM64_JMG(k,:));
hold on
legent{k}=[num2str(itlim(k))];
end
hold off
legend(legent)
xlabel('SNR (dB)')
ylabel('BER')
title({'Comparison of BER for LDPC for different iteration limits','QAM64'})

rmpath(genpath(pwd))