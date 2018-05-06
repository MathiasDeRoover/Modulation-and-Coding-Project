function [outStream] = TruncateAndSample(inStream,ftaps,periodToSample,streamFrequency,SCO,timeShift)
%TRUNCATEANDSAMPLE Truncate and sample the inStream
%   Throw away 2*ftaps asmples on both sides of signal and sample it at
%   sample frequency 1/periodToSample.
truncStream     = inStream(2*ftaps+1:end-2*ftaps);
truncStream2    = inStream(2*ftaps+1:end);
t               = ((1:numel(truncStream))-1)/streamFrequency;
t2              = ((1:numel(truncStream2))-1)/streamFrequency;
tq              = ((1:(t(end)/periodToSample)+1)-1)*(periodToSample * (1 + SCO)) + timeShift;
outStream       = interp1(t2,truncStream2,tq)';
end

