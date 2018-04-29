%% Initializations
clear
close all
clc
addpath(genpath(pwd))

%% Generate H matrix
c_nodes     = 128;                  % Block length
v_nodes     = 256;                  % Coded block length
SNR         = linspace(-20,20,200);  % SNR in dB for IdealChannel_exec
% SNR         = db2mag(SNRdb);
H           = makeLdpc(c_nodes,v_nodes,0,1,3);  % Create initial parity check matrix

BPSK_JMG        = zeros(size(SNR));
BPSK_Guylian    = zeros(size(SNR));
UncodedBER      = zeros(size(SNR));


%% Generate bitstream
N           = 768*20;   % Number of symbols
bps         = 1;    % Bits per symbol

wait_bar = waitbar(0,'Simulating');
itlim = [1 10 50 100];
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

    %% Send through channel
%     if j==1
%         receivedUncodedbpsk     = IdealChannel_exec(bitStream,SNR(i),'BPSK','det',bitStream);
%     end
    %     receivedUncodedqpsk     = IdealChannel_exec(bitStream,SNR(i),'QPSK','det');
    %     receivedUncoded16qam     = IdealChannel_exec(bitStream,SNR(i),'16QAM','det');
    %     receivedUncoded64qam     = IdealChannel_exec(bitStream,SNR(i),'64QAM','det');
    %     
    %     receivedCodedbpsk     = IdealChannel_exec(bitSCoded,SNR(i),'BPSK','det');
    %     receivedCodedqpsk     = IdealChannel_exec(bitSCoded,SNR(i),'QPSK','det');
    %     receivedCoded16qam     = IdealChannel_exec(bitSCoded,SNR(i),'16QAM','det');
    %     receivedCoded64qam     = IdealChannel_exec(bitSCoded,SNR(i),'64QAM','det');

        receivedStream1 = IdealChannel_exec(bitSCoded,SNR(i),'BPSK','det',bitStream);
        receivedStream2 = IdealChannel_exec(bitSCoded,SNR(i),'QPSK','det',bitStream);
        receivedStream3 = IdealChannel_exec(bitSCoded,SNR(i),'16QAM','det',bitStream);
        receivedStream4 = IdealChannel_exec(bitSCoded,SNR(i),'64QAM','det',bitStream);

    %% Decode bitstream

    %     hardDecodedJMG      = LDPC_decoder_hard_biased( receivedStream, Hs, 10 );
%         if j==1
%         hardDecodedGuylian      = hardDecoderLDPC_G2( receivedStream, Hs, c_nodes, v_nodes );
%         end
        hardDecoded1    = LDPC_decoder_hard_biased( receivedStream1, Hs, itlim(j) );
        hardDecoded2    = LDPC_decoder_hard_biased( receivedStream2, Hs, itlim(j) );
        hardDecoded3    = LDPC_decoder_hard_biased( receivedStream3, Hs, itlim(j) );
        hardDecoded4    = LDPC_decoder_hard_biased( receivedStream4, Hs, itlim(j) );
        
    %     bitRecoveredHardbpsk    = LDPC_decoder_hard( receivedCodedbpsk, Hs, 10 );
    % %     bitRecoveredHardbpsk    = LDPC_decoder_hard_lite( receivedCodedbpsk, Hs);
    % 
    %     bitRecoveredHardqpsk    = LDPC_decoder_hard( receivedCodedqpsk, Hs, 10 );
    % %     bitRecoveredHardqpsk    = LDPC_decoder_hard_lite( receivedCodedqpsk, Hs);
    % 
    %     bitRecoveredHard16qam    = LDPC_decoder_hard( receivedCoded16qam, Hs, 10 );
    % %     bitRecoveredHard16qam    = LDPC_decoder_hard_lite( receivedCoded16qam, Hs);
    %     
    %     bitRecoveredHard64qam    = LDPC_decoder_hard( receivedCoded64qam, Hs, 10 );
    % %     bitRecoveredHard64qam    = LDPC_decoder_hard_lite( receivedCoded64qam, Hs);

    %% Calculate bit error
        [~,BPSK_JMG(j,i)]     = biterr(hardDecoded1,bitStream);
        [~,QPSK_JMG(j,i)]     = biterr(hardDecoded2,bitStream);
        [~,QAM16_JMG(j,i)]     = biterr(hardDecoded3,bitStream);
        [~,QAM64_JMG(j,i)]     = biterr(hardDecoded4,bitStream);
        
        
        if j==1
%         [~,BPSK_Guylian(i)] = biterr(hardDecodedGuylian,bitStream);
%         [~,UncodedBER(i)]   = biterr(receivedUncodedbpsk,bitStream);
        end

    %     [~,BPSK_BERH(i)] = biterr(bitRecoveredHardbpsk,bitStream);
    %     [~,BPSK_BERN(i)] = biterr(receivedUncodedbpsk,bitStream);
    %     
    %     [~,QPSK_BERH(i)] = biterr(bitRecoveredHardqpsk,bitStream);
    %     [~,QPSK_BERN(i)] = biterr(receivedUncodedqpsk,bitStream);
    %     
    %     [~,QAM16_BERH(i)] = biterr(bitRecoveredHard16qam,bitStream);
    %     [~,QAM16_BERN(i)] = biterr(receivedUncoded16qam,bitStream);
    %     
    %     [~,QAM64_BERH(i)] = biterr(bitRecoveredHard64qam,bitStream);
    %     [~,QAM64_BERN(i)] = biterr(receivedUncoded64qam,bitStream);
        waitbar(((j-1)*numel(SNR)+i)/numel(SNR)/numel(itlim),wait_bar);
    end
end
close(wait_bar);

%% Plots
% figure
% chosenlim=10;
% semilogy(SNR, BPSK_JMG(find(itlim==chosenlim),:));
% hold on
% semilogy(SNR, BPSK_Guylian)
% semilogy(SNR, UncodedBER)
% hold off
% legend('Our BPSK','Guylian BPSK','Uncoded BER')
% xlabel('SNR (dB)')
% ylabel('BER')
% title(['Comparison of BER for LDPC with iteration limit = ' num2str(chosenlimit)])



% BPSK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
legent=[] ;
for k=1:numel(itlim)
semilogy(SNR, BPSK_JMG(k,:));
legent{k}=[num2str(itlim(k))];
end
hold off
legend(legent)
xlabel('SNR (dB)')
ylabel('BER')
title({'Comparison of BER for LDPC for different iteration limits','BPSK'})

% QPSK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
legent=[] ;
for k=1:numel(itlim)
semilogy(SNR, QPSK_JMG(k,:));
legent{k}=[num2str(itlim(k))];
end
hold off
legend(legent)
xlabel('SNR (dB)')
ylabel('BER')
title({'Comparison of BER for LDPC for different iteration limits','QPSK'})

% 16QAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
legent=[] ;
for k=1:numel(itlim)
semilogy(SNR, QAM16_JMG(k,:));
legent{k}=[num2str(itlim(k))];
end
hold off
legend(legent)
xlabel('SNR (dB)')
ylabel('BER')
title({'Comparison of BER for LDPC for different iteration limits','QAM16'})

% 64QAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
legent=[] ;
for k=1:numel(itlim)
semilogy(SNR, QAM64_JMG(k,:));
legent{k}=[num2str(itlim(k))];
end
hold off
legend(legent)
xlabel('SNR (dB)')
ylabel('BER')
title({'Comparison of BER for LDPC for different iteration limits','QAM64'})

% figure
% subplot(2,2,1)
% semilogy(SNR, BPSK_BERH)
% hold on
% semilogy(SNR, BPSK_BERN)
% hold off
% legend('Not coded','Hard')
% xlabel('SNR (dB)')
% ylabel('BER')
% title('BPSK')
% 
% subplot(2,2,2)
% semilogy(SNR, QPSK_BERH)
% hold on
% semilogy(SNR, QPSK_BERN)
% hold off
% legend('Not coded','Hard')
% xlabel('SNR (dB)')
% ylabel('BER')
% title('QPSK')
% 
% subplot(2,2,3)
% semilogy(SNR, QAM16_BERH)
% hold on
% semilogy(SNR, QAM16_BERN)
% hold off
% legend('Not coded','Hard')
% xlabel('SNR (dB)')
% ylabel('BER')
% title('QAM16')
% 
% subplot(2,2,4)
% semilogy(SNR, QAM64_BERH)
% hold on
% semilogy(SNR, QAM64_BERN)
% hold off
% legend('Not coded','Hard')
% xlabel('SNR (dB)')
% ylabel('BER')
% title('QAM64')


rmpath(genpath(pwd))