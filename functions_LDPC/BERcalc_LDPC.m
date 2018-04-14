function [ BERratio ] = BERcalc_LDPC(N,ftaps,modu,M,fs,g,EbN0_array,current_loop,number_of_loops )
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
if mod(numel(bitStream),bps) ~= 0
    error('Number of bits not a multiple of bps');
end

symStream = mapping(bitStream, bps, modulation);
supStream = upsample(symStream,M);
sgStream  = conv(supStream,g);



%% Ideal Channel
SignalEnergy = (trapz(abs(sgStream(ftaps+1:end-(ftaps-1))).^2))*(1/fs);        % Total signal energy (this is given in slides) %Trapz is integral approx.
Eb = SignalEnergy /(numel(bitStream));
Eb = Eb/2;

%%% Calculating BER for every EbN0 in EbN0_array
BERratio = zeros(1,numel(EbN0_array));
wait_bar = waitbar(0,['Executing loop ',num2str(current_loop),'/',num2str(number_of_loops),' ',modu,'...']);
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
    
    [~,BERratio(i)]=biterr(recStream,bitStream);
    waitbar(i/numel(EbN0_array),wait_bar);
end
close(wait_bar)
end

