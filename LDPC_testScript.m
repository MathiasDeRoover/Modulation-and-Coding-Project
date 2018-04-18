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
% % noise = rand(1,numel(bitStream_enc))'>0.99;
% 
% % bitStream_enc'
% 
% % bitStream_enc = mod(bitStream_enc+noise,2);
% bitStream_enc(5:10:end) = ~bitStream_enc(5:10:end);         % Add some biterrors manually
% 
% tic
% % bitStream_rec = LDPC_decoder_hard( bitStream_enc, newH );
% toc
% 
% tic
% bitStream_rec = LDPC_decoder_hard_lite( bitStream_enc, H );
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
N                           = 100;
c_length                    = 5;
v_length                    = 10;
bitStream                   = CreateBitStream(N,c_length);
%bitStream = [1 0 0 1 1]';
H0                          = makeLdpc(c_length,v_length,0,1,3);
[bitStream_enc,newH]        = LDPC_encoder_lite( bitStream, H0 ); %Was prev H

bitStream_chan = real(IdealChannel_exec(bitStream_enc,200,'BPSK','no_demodu'));
figure
histogram(abs(bitStream_chan),0:0.05:2)

tic
bitStream_rec = LDPC_decoder_soft_log_BPSK( bitStream_chan,H,std(bitStream_chan(bitStream_chan>0)) );
toc

%% Plotting results

figure
stem(bitStream ~= bitStream_rec);

rmpath(genpath(pwd))