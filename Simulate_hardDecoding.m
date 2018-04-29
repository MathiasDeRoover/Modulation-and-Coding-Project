%% Initializations
clear
close all
clc
addpath(genpath(pwd))

%% Run Simulations
N                           = 100;   % Number of Points
c_length                    = 128;
v_length                    = 256;
bitStream                   = CreateBitStream(N,c_length);
H0                          = makeLdpc(c_length,v_length,0,1,3);

bitStream_blk               = reshape(bitStream,c_length,[]);
[bitStream_cod_blk,newH]    = makeParityChk(bitStream_blk, H0, 0);  % Create parity check bits and reshape H
bitStream_cod_blk           = [bitStream_cod_blk;bitStream_blk];    % Unite parity check bits and message
bitStream_cod               = reshape(bitStream_cod_blk,[],1);

modu = {'BPSK','QPSK','16QAM','64QAM'};
f1 = figure;
f2 = figure;

SNR = linspace(-20,20,100);
for m = 1:4
    noDecBer = zeros(size(SNR));
    hard1ber = zeros(size(SNR));
    hard1t   = zeros(size(SNR));
    
    wait_bar = waitbar(0,'please wait...');
    for i = 1:numel(SNR)
        
        bitStream_chan = real(IdealChannel_exec(bitStream_cod,SNR(i),modu{m},'det'));
        
        bitStream_chan_det_block    = reshape(bitStream_chan,v_length,[]);
        bitStream_rec_block         = bitStream_chan_det_block(end-c_length+1:end,:);
        bitStream_rec               = reshape(bitStream_rec_block,[],1);
        [~,noDecBer(i)] = biterr(bitStream_rec,bitStream);
        
        tic
        bitStream_rec   = LDPC_decoHardVec( bitStream_chan, newH ,10);
        hard1t(i) = toc ;
        [~,hard1ber(i)] = biterr(bitStream_rec,bitStream);
        
        waitbar(i/numel(SNR),wait_bar);
    end
    close(wait_bar)
    
    %% Plotting results
    figure(f1)
    subplot(2,2,m)
    semilogy(SNR,noDecBer)
    hold on
    semilogy(SNR,hard1ber)
    hold off
    ylabel('BER')
    xlabel('SNR')
    legend('No decoding','hard decode')
    title(modu(m))
    
    figure(f2)
    subplot(2,2,m)
    hold on
    plot(SNR,hard1t)
    hold off
    ylabel('timing [s]')
    xlabel('SNR')
    legend('hard decode')
    title(modu{m})
end

rmpath(genpath(pwd))



% figure
% semilogy(SNRdb,BPSKarray{1},'b');
% hold on
% plot([4 6 8],[2e-2 3e-3 3e-4]./log2(2),'bo');
% semilogy(SNRdb,QPSKarray{1},'Color',[0.80 0.40 0.0]);
% plot([1 5 7],[1e-1 1e-2 1e-3]./log2(4),'o','Color',[0.80 0.40 0.0]);
% semilogy(SNRdb,Qam16array{1},'r');
% plot([9 10 12],[2e-2 5e-3 1e-3]./log2(16),'ro');
% semilogy(SNRdb,Qam64array{1},'k');
% plot([12 14 18],[5e-2 1e-2 1e-4]./log2(64),'ko');
% hold off
% xlabel('SNR=Eb/N0(dB)')
% ylabel('BER')
% legend('BPSK','BPSKth','QPSK','QPSKth','Qam16','Qam16th','Qam64','Qam64th')
% title(['BER for different constellations, upsampling rate = ' num2str(M(1))])