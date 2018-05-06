function [outStream] = Gardner(inStream, K, windStream, ftaps, windFrequency, periodToSample,SCO, timeShift)
%GARDNER Summary of this function goes here
%   Detailed explanation goes here
N = numel(inStream)/2;
e = zeros(N,1);
windStream  = windStream(2*ftaps+1:end);
windTime    = ((1:numel(windStream))-1)/windFrequency;
epsOld  = 0;
eps     = 0;
outStream = zeros(N,1);
for n = 1:N-1
    timeSample1     = (n-1)*((periodToSample*(1 + SCO)) + timeShift - epsOld);
    timeHalfSample  = (n-1/2)*((periodToSample*(1 + SCO)) + timeShift - eps);
    timeSample2     = (n)*((periodToSample*(1 + SCO)) + timeShift - eps);
   
    sample1         = interp1(windTime,windStream,timeSample1);
    halfSample      = interp1(windTime,windStream,timeHalfSample);
    sample2         = interp1(windTime,windStream,timeSample2);
    
    e(n+1)          = e(n) + 2*K*real( halfSample * (conj(sample2) - conj(sample1)) );
    eps             = e(n+1);
    epsOld          = e(n);
    outStream(n)    = sample1;
end
outStream(N)        = sample2;
end

