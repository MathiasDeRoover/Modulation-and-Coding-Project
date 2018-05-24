function [outStream] = TruncateAndSample(inStream,ftaps,periodToSample,streamFrequency,SCO,timeShift)
%TRUNCATEANDSAMPLE Truncate and sample the inStream
% DISCLAIMER: extra lines in comments have been added so Johan can
% understand the code. These can be removed in the final version.
%   Throw away 2*ftaps samples on both sides of signal and sample it at
%   sample frequency 1/periodToSample.
% truncStream = inStream; %TEST
truncStream    = inStream(2*ftaps+1:end);                           % Last points are not removed to prevent interpolation errors. 
tend            = (numel(inStream)-1-4*ftaps)/streamFrequency;
% tend            = (numel(truncStream)-1-2*ftaps)/streamFrequency; % Last sample at 2*ftaps before end in case of truncStream
t2              = ((1:numel(truncStream))-1)/streamFrequency.';
tq              = ((1:(tend/periodToSample)+1)-1)*(periodToSample * (1 + SCO)) + timeShift;
% tq            = (0:periodToSample:tend) *(1 + SCO) + timeShift; % Analoog
outStream       = interp1(t2,truncStream,tq)';
outStream = conj(outStream);

% truncStream2 = inStream(2*ftaps+1:end-2*ftaps);
% t3              = ((1:numel(truncStream2))-1)/streamFrequency.';
% figure
% plot(t3,real(truncStream2),'*r',tq,real(outStream),'-r');
% hold on
% plot(t3,imag(truncStream2),'*b',tq,-imag(outStream),'-b');
% legend('Real input','Real interp1','Imag input','Imag interp1')

%interp1:truncStream=F(t2)-> find outstream=F(tq)
end

