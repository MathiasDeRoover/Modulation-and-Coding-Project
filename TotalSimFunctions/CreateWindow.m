function [g,g_min] = CreateWindow(symbolPeriod, sampleFrequency, ftaps, beta)
%CREATEWINDOW Summary of this function goes here
%   Detailed explanation goes here 
H = RRCFilter(symbolPeriod,sampleFrequency,beta,ftaps )'; % Creating the window
h = ifft(H,'symmetric');                                  % Transforming to the time domain
h = h/h(1);                                               % Normalizing the window in the time domain
H = fft(h);                                               % Transforming the normalized window back to the frequency domain

G = sqrt(H);                                              % G is the square root of H so that G*G = H (after the convolutions)
g = fftshift(ifft(G,'symmetric'));                        % Transform G to the time domain. Fftshift is needed to get a proper raised cosine
g_min=fliplr(g);
end

