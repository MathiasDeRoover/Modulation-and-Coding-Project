function [ ] = plotFilter( fs,h,H,G,M,ftaps )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%% Plotting the window
plot_f = -fs/2:fs/2/ftaps:fs/2-fs/2/ftaps;                  % Frequency axis window 
plot_t = (-floor(numel(h)/2):1:floor(numel(h)/2)-1)*1/fs;   % Time axis window
plot_H = fftshift(real(H));                                 % H ready to plot
plot_h = fftshift(ifft(H,'symmetric'));                     % h ready to plot
plot_g = fftshift(ifft(G,'symmetric'));                     % g ready to plot
plot_gmin = fliplr(plot_g);                                 % g(-t) ready to plot

figure()

subplot(2,2,1) %plot1
plot(plot_f,plot_H);
title('H in f domain')
xlabel('Frequency [Hz]')
ylabel('H(f)')

subplot(2,2,2) %plot2
hold on
plot(plot_t,plot_h);
stem(plot_t(1:M:end),plot_h(1:M:end))
title('H in t domain')
xlabel('time [s]')
ylabel('h(t)')
hold off

subplot(2,2,3) %plot3
plot(plot_t,plot_g);
title('g(t)')
xlabel('time [s]')
ylabel('g(t)')

subplot(2,2,4) %plot4
plot(plot_t,plot_gmin);
title('g(-t)')
xlabel('time [s]')
ylabel('g(-t)')

end

