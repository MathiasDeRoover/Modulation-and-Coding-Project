function [] = Validation( shsStream,symStream,bitStream,recStream,T )
%Function to validate (plot) channel performance
%%% Plotting the received signals
figure
hold on
plot(real(shsStream),imag(shsStream),'*')
plot(real(symStream),imag(symStream),'*')
hold off
title('Received symbols VS Send symbols')
%%%

%%% Plotting the end result (Correct if it is zero everywhere)
plot_t2 = (0:numel(recStream)-1)*T;
figure()
hold on
stem(plot_t2,bitStream-recStream);
hold off
title('Original bitsream minus Recovered bitstream')
xlabel('time [s]')
ylabel('Correct if = 0')
%%%

end

