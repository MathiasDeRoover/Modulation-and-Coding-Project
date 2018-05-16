function [H] = RRCFilter( T,fs,beta,ftaps )
%RRCFilter Calculate the RCC filter
%   Calculate the RCC filter frequency spectrum in such a way that is is
%   the same length as the upsampled symbolstream frequencySpectrum. Needs
%   the symbol period, sample frequency, roll off factor and amount of
%   causal filter taps.

x1 = (1-beta)/(2*T);
x2 = (1+beta)/(2*T);
f = linspace(-fs/2,fs/2,2*ftaps+1);

H =     T                                       .* ((0<=abs(f)) & (abs(f)<x1)) +...
        T/2.*(1+cos(pi*T/beta.*(abs(f)-x1)))    .* ((x1<=abs(f)) & (abs(f)<=x2)) +...
        0                                       .* (abs(f)>x2);
H = ifftshift(H);
H = H(1:end-1);      % Drop the last sample so that fs/2 <= f < fs/2
end

