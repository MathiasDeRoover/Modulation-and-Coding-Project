function [outStream] = TruncateAndSample(inStream,ftaps,periodToSample,streamFrequency,SCO,timeShift)
%TRUNCATEANDSAMPLE Truncate and sample the inStream
% DISCLAIMER: extra lines in comments have been added so Johan can
% understand the code. These can be removed in the final version.
%   Throw away 2*ftaps samples on both sides of signal and sample it at
%   sample frequency 1/periodToSample.
truncStream    = inStream(2*ftaps+1:end);                           % Last points are not removed to prevent interpolation errors. 
tend            = (numel(inStream)-1-4*ftaps)/streamFrequency;
% tend            = (numel(truncStream)-1-2*ftaps)/streamFrequency; % Last sample at 2*ftaps before end in case of truncStream
t2              = ((1:numel(truncStream))-1)/streamFrequency;
tq              = ((1:(tend/periodToSample)+1)-1)*(periodToSample * (1 + SCO)) + timeShift;
% tq            = (0:periodToSample:tend) *(1 + SCO) + timeShift; % Analoog
outStream       = interp1(t2,truncStream,tq)';
end

