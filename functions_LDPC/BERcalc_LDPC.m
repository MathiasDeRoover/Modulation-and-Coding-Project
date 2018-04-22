function [ BERratio ] = BERcalc_LDPC(N,ftaps,modu,M,fs,g,EbN0_array,iteration_limit )
%Based on modulation type and SNR ratio, the function returns the BER.

%% Tranceiver
switch modu
    case 'BPSK'
        modulation  = 'psk';
        bps         = 1;
    case 'QPSK'
        modulation  = 'psk';
        bps         = 2;
    case '16QAM'
        modulation  = 'qam';
        bps         = 4;
    case '64QAM'
        modulation  = 'qam';
        bps         = 6;
    otherwise
        error([modu,' is an unsupported constellation'])
end

bitStream = CreateBitStream(N,bps);
bitStreamCheck = bitStream;

if mod(numel(bitStream),bps) ~= 0
    error('Number of bits not a multiple of bps');
end

%Encoding
% v_length                    = N;
% c_length                    = N/2;
v_length                    = 256;
c_length                    = 128;

zerosToPad = c_length - mod(length(bitStream),c_length);
if (zerosToPad ~= 0)
%     disp('We need to pad some zeros')
    bitStreamPad = [bitStream; zeros(zerosToPad,1)];
    while mod(numel(bitStream),bps) ~= 0                % Infinite loops should not happen, but be careful!
        bitStreamPad = [bitStreamPad; zeros(c_length,1)];
        zerosToPad = zerosToPad + c_length;
    end
end

H0                          = makeLdpc(c_length,v_length,0,1,3);            % Create initial parity check matrix of size 128 x 256

tic
bitStream_blk               = reshape(bitStreamPad,c_length,[]);
[bitStream_cod_blk,newH]    = makeParityChk(bitStream_blk, H0, 0);          % Create parity check bits and reshape H
bitStream_cod_blk           = [bitStream_cod_blk;bitStream_blk];            % Unite parity check bits and message
bitStream_cod               = reshape(bitStream_cod_blk,[],1);      

%Rest of channel
% Remove padding after encoding
% bitStream_cod = bitStream_cod(1:end-zerosToPad);

symStream = mapping(bitStream_cod, bps, modulation);
supStream = upsample(symStream,M);
sgStream  = conv(supStream,g);



%% Ideal Channel
%Do we have to use sgStream or bitStream? (encoded vs unencoded)
%SignalEnergy = (trapz(abs(sgStream(ftaps+1:end-(ftaps-1))).^2))*(1/fs); % Total signal energy (this is given in slides) %Trapz is integral approx.
SignalEnergy = (trapz(abs(bitStream).^2))*(1/fs); 
Eb = SignalEnergy /(numel(bitStream));
Eb = Eb/2;

%%% Calculating BER for every EbN0 in EbN0_array
BERratio = zeros(1,numel(EbN0_array));
for i = 1:numel(EbN0_array)
    EbN0 = EbN0_array(i);
    
    N0 = Eb/EbN0;
    NoisePower = 2*N0*fs;
    noise = sqrt(NoisePower/2) * (randn(numel(sgStream),1) + 1i*randn(numel(sgStream),1));
    sgStream_n = sgStream + noise;
    
    
    %% Receiver
    gmin=fliplr(g);                                         % Converting g(t) to g(-t) to get matched filter
    switch modulation                                       % Windowing downStream in the frequencyDomain with g(-t)
        case 'pam'
            sggStream = real(conv(sgStream_n,gmin));        % Taking real because noise is complex and BPSK signal is real.
        otherwise
            sggStream = conv(sgStream_n,gmin);              % Noise and qam signal are complex.
    end
    
    sggStream = sggStream(2*ftaps+1 : end-2*(ftaps-1));     % Dropping access data that originates from convolutions
    shsStream = sggStream(1:M:end);                         % Sampling at nT
    recStream = demapping(shsStream, bps, modulation);      % Demapping
    
    %Decoder
    bitStream_rec  = LDPC_decoder_hard(recStream, newH ,iteration_limit);
    
    % Receiver knows how much padding is added
    bitStream_rec   = bitStream_rec(1:end-zerosToPad);

%     result = bitStream_rec ~= bitStreamCheck;
%     figure % MAKE SURE YOU COMPUTE FOR 1 SNR
%     stem(result)

    %BER
    [~,BERratio(i)]=biterr(bitStream_rec,bitStreamCheck);

end

end

