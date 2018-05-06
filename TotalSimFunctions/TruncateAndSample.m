function [outStream] = TruncateAndSample(inStream,ftaps,periodToSample,streamFrequency,SCO,timeShift)
%TRUNCATEANDSAMPLE Truncate and sample the inStream
%   Throw away 2*ftaps asmples on both sides of signal and sample it at
%   sample frequency 1/periodToSample.
truncStream    = inStream(2*ftaps+1:end);
tend            = (numel(inStream)-1-4*ftaps)/streamFrequency;
t2              = ((1:numel(truncStream))-1)/streamFrequency;
tq              = ((1:(tend/periodToSample)+1)-1)*(periodToSample * (1 + SCO)) + timeShift;
outStream       = interp1(t2,truncStream,tq)';
end

