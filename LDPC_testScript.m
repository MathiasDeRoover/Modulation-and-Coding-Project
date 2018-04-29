%% Initializations
clear
close all
clc
addpath(genpath(pwd))

H = [   1 1 0 1 1 0 0 1 0 0;        % Define Parity Check matrix H:
    0 1 1 0 1 1 1 0 0 0;        %   This matrix is sparse,
    0 0 0 1 0 0 0 1 1 1;        %   it exists out of a lot
    1 1 0 0 0 1 1 0 1 0;        %   of zeros and a few ones.
    0 0 1 0 0 1 0 1 0 1;   ];   %

%% Step 1
% N = 1e2;
% bps = 4;
% bitStream = CreateBitStream(N,bps);
% [bitStream_enc,newH] = LDPC_encoder_lite( bitStream, H );
% noise = rand(1,numel(bitStream_enc))'>0.99;

% bitStream_enc'

% bitStream_enc = mod(bitStream_enc+noise,2);
% bitStream_enc(5:10:end) = ~bitStream_enc(5:10:end);         % Add some biterrors manually

%tic
% bitStream_rec = LDPC_decoder_hard( bitStream_enc, newH );
%toc

% tic
% bitStream_rec = LDPC_decoder_hard_lite( bitStream_enc, H );
% toc

% tic
% bitStream_rec = block_decoder( bitStream_enc, H );
% toc

%% Step 2
% N                           = 50;
% c_length                    = 128;
% v_length                    = 256;
% bitStream                   = CreateBitStream(N,c_length);
% H0                          = makeLdpc(c_length,v_length,0,1,3);            % Create initial parity check matrix of size 128 x 256
%
% tic
% bitStream_blk               = reshape(bitStream,c_length,[]);
% [bitStream_cod_blk,newH]    = makeParityChk(bitStream_blk, H0, 0);          % Create parity check bits and reshape H
% bitStream_cod_blk           = [bitStream_cod_blk;bitStream_blk];            % Unite parity check bits and message
% bitStream_cod               = reshape(bitStream_cod_blk,[],1);
% toc
%
% tic
% bitStream_rec               = LDPC_decoder_hard( bitStream_cod, newH,10 );     % Decode
% toc

%% Step 3


% IMPORTANT: this code is incorrect because it uses the LDPC_encoder_lite
% in combination with a 128x256 H matrix. This will result in a H matrix
% with non integer values. The code has been corrected in branch
% Fixing_softDecoding




N                           = 60;
c_length                    = 128;
v_length                    = 256;
bitStream                   = CreateBitStream(N,c_length);
%bitStream = [1 0 0 1 1]';
H0                          = makeLdpc(c_length,v_length,0,1,3);

bitStream_blk               = reshape(bitStream,c_length,[]);
[bitStream_cod_blk,newH]    = makeParityChk(bitStream_blk, H0, 0);          % Create parity check bits and reshape H
bitStream_cod_blk           = [bitStream_cod_blk;bitStream_blk];            % Unite parity check bits and message
bitStream_cod               = reshape(bitStream_cod_blk,[],1);

modu = {'BPSK','QPSK','16QAM','64QAM'};
f1 = figure;
f2 = figure;

for m = 1:4
    SNR      = linspace(-20,20,100);
    noDecBer    = zeros(size(SNR));
    hard2ber    = zeros(size(SNR));
    hard2t      = zeros(size(SNR));
    hard1ber    = zeros(size(SNR));
    hard1t      = zeros(size(SNR));
    
    wait_bar = waitbar(0,'please wait...');
    for i = 1:numel(SNR)
        
        bitStream_chan = real(IdealChannel_exec(bitStream_cod,SNR(i),modu{m},'det'));
        
        bitStream_chan_det_block    = reshape(bitStream_chan,v_length,[]);
        bitStream_rec_block         = bitStream_chan_det_block(end-c_length+1:end,:);
        bitStream_rec               = reshape(bitStream_rec_block,[],1);
        [~,noDecBer(i)] = biterr(bitStream_rec,bitStream);
        
        tic
        bitStream_rec   = LDPC_decoder_hard_biased( bitStream_chan, newH, 10);
        hard2t(i) = toc;
        [~,hard2ber(i)] = biterr(bitStream_rec,bitStream);
        
        tic
        bitStream_rec   = LDPC_decoder_hard( bitStream_chan, newH ,10);
        hard1t(i) = toc ;
        [~,hard1ber(i)] = biterr(bitStream_rec,bitStream);
        
        waitbar(i/numel(SNR),wait_bar);
    end
    close(wait_bar)
    % tic
    % bitStream_rec = LDPC_decoder_soft_BPSK( bitStream_chan,newH,std(bitStream_chan(bitStream_chan>0)) );
    % toc
    
    %% Plotting results
    figure(f1)
    subplot(2,2,m)
    semilogy(SNR,noDecBer)
    hold on
    semilogy(SNR,hard1ber)
    semilogy(SNR,hard2ber)
    hold off
    ylabel('BER')
    xlabel('Max iterations at SNR = 5')
    legend('No decoding','hard decode','hard decode biased')
    title(modu(m))
    
    figure(f2)
    subplot(2,2,m)
    hold on
    plot(SNR,hard1t)
    plot(SNR,hard2t)
    hold off
    ylabel('timing [s]')
    xlabel('Max iterations at SNR = 5')
    legend('hard decode','hard decode biased')
    title(modu{m})
end

rmpath(genpath(pwd))