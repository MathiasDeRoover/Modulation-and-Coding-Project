%% Initializations
clear
close all
clc
addpath(genpath(pwd))

H = [   1 1 0 1 1 0 0 1 0 0;        % Define Parity Check matrix H:
        0 1 1 0 1 1 1 0 0 0;        %   This matrix is sparse,
        0 0 0 1 0 0 0 1 1 1;        %   it exists out of a lot
        1 1 0 0 0 1 1 0 1 0;        %   of zeros and ones.
        0 0 1 0 0 1 0 1 0 1;   ];   %

%% Execution

bitStream = CreateBitStream(1,5);
bitStream_enc = LDPC_encoder_lite( bitStream, H );
%noise = rand(1,numel(bitStream_enc))'>0.8;

bitStream_enc'

%bitStream_enc = mod(bitStream_enc+noise,2);
bitStream_enc(3:10:end) = ~bitStream_enc(3:10:end);
bitStream_rec = LDPC_decoder_hard_lite( bitStream_enc, H );

%% Plotting results

figure
stem(bitStream ~= bitStream_rec);

rmpath(genpath(pwd))