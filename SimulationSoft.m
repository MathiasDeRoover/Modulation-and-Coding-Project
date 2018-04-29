%% Initializations
clear
close all
clc
addpath(genpath(pwd))

%% Create H
N                           = 60;
c_length                    = 128;
v_length                    = 256;
bitStream                   = CreateBitStream(N,c_length);
H0                          = makeLdpc(c_length,v_length,0,1,3);

% c_length                    = 5;
% v_length                    = 10;
% H0 = [   1 1 0 1 1 0 0 1 0 0;        % Define Parity Check matrix H:
%     0 1 1 0 1 1 1 0 0 0;        %   This matrix is sparse,
%     0 0 0 1 0 0 0 1 1 1;        %   it exists out of a lot
%     1 1 0 0 0 1 1 0 1 0;        %   of zeros and a few ones.
%     0 0 1 0 0 1 0 1 0 1;   ];   %
% bitStream                   = CreateBitStream(N,c_length);

%% Encode
bitStream_blk               = reshape(bitStream,c_length,[]);
[bitStream_cod_blk,newH]    = makeParityChk(bitStream_blk, H0, 0);          % Create parity check bits and reshape H
bitStream_cod_blk           = [bitStream_cod_blk;bitStream_blk];            % Unite parity check bits and message
bitStream_cod               = reshape(bitStream_cod_blk,[],1);

f1 = figure;
f2 = figure;

SNR        = linspace(-20,20,100);
noDecBer   = zeros(size(SNR));
hardBer    = zeros(size(SNR));
hardT      = zeros(size(SNR));
softBer    = zeros(size(SNR));
softT      = zeros(size(SNR));
softLogBer = zeros(size(SNR));
softLogT   = zeros(size(SNR));

wait_bar = waitbar(0,'please wait...');
for i = 1:numel(SNR)
    
    bitStream_chan = real(IdealChannel_exec(bitStream_cod,SNR(i),'BPSK','no_det'));
    sigma = std(real(bitStream_chan(real(bitStream_chan)>0)));
    bitStream_chan_det = real(bitStream_chan) > 0.5;
    
    bitStream_chan_det_block    = reshape(bitStream_chan_det,v_length,[]);
    bitStream_rec_block         = bitStream_chan_det_block(end-c_length+1:end,:);
    bitStream_rec               = reshape(bitStream_rec_block,[],1);
    [~,noDecBer(i)] = biterr(bitStream_rec,bitStream);
    
    tic
    bitStream_rec   = LDPC_decoder_hard( bitStream_chan_det, newH ,200);
    hardT(i) = toc ;
    [~,hardBer(i)] = biterr(bitStream_rec,bitStream);
    
    tic
    bitStream_rec   = LDPC_decoder_soft2(real(bitStream_chan), newH, sigma, 200);
    softT(i) = toc ;
    [~,softBer(i)] = biterr(bitStream_rec,bitStream);
    
    tic
    bitStream_rec   = LDPC_decoder_softLog2(real(bitStream_chan), newH, sigma, 2);
%     bitStream_rec   = hardDecoderLDPC_G2( bitStream_chan_det, newH, 128, 256 );
    softLogT(i) = toc ;
    [~,softLogBer(i)] = biterr(bitStream_rec,bitStream);
    
    waitbar(i/numel(SNR),wait_bar);
end
close(wait_bar)
% tic
% bitStream_rec = LDPC_decoder_soft_BPSK( bitStream_chan,newH,std(bitStream_chan(bitStream_chan>0)) );
% toc

%% Plotting results
figure(f1)
semilogy(SNR,noDecBer)
hold on
semilogy(SNR,hardBer)
semilogy(SNR,softBer)
semilogy(SNR,softLogBer)
hold off
ylabel('BER')
xlabel('Max iterations at SNR = 5')
legend('No decoding','hard','soft','log')
title('BPSK')

figure(f2)
hold on
plot(SNR,hardT)
plot(SNR,softT)
plot(SNR,softLogT)
hold off
ylabel('timing [s]')
xlabel('Max iterations at SNR = 5')
legend('hard','soft','log')
title('BPSK')

rmpath(genpath(pwd))