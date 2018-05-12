function [modulation,bps] = ModuToModulation(modu)
%MODUTOMODULATION Transform modu into modulation and bps
%   Transform modu into modulation and bps
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
end

