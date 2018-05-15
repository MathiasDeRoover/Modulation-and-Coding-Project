function [noisePower] = CalcNoisePower(bitStream, symStream,streamT,SNRdB,fs,ftaps)
%CALCNOISEPOWER Calculate the noisepower based on stream and SNR
%   Calculate the noisepower based on stream and SNR
SignalEnergy = (trapz(abs(symStream(ftaps+1:end-(ftaps-1))).^2))*streamT;
Eb = SignalEnergy / numel(bitStream);                   % Energy for a single bit? see slides
Eb = Eb/2;                                           % The power of a bandpass signal is equal to the power of its complex envelope divided by 2
EbN0 = db2mag(SNRdB*2);                              % Variable SNR; factor 2 because we are working with powers, not amplitudes: powerdB = 10*log10(power) vs amplitudedB = 20*log10(amplitude)
N0 = Eb/EbN0;
noisePower = 2*N0*fs;
end

