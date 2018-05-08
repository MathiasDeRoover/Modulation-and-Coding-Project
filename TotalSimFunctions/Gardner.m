function [outStream] = Gardner(inStream, K, windStream, ftaps, windFrequency, periodToSample,SCO, timeShift)
%GARDNER Summary of this function goes here
%   Detailed explanation goes here
timeShift = 0;              % You are not supposed to know this here I think
SCO = 0;                    % You are not supposed to know this here I think
N = numel(inStream)/2;
e = zeros(N,1);
windStream  = windStream(2*ftaps+1:end);
windTime    = ((1:numel(windStream))-1)/windFrequency;
epsOld  = 0;
eps     = 0;
outStream = zeros(N,1);
for n = 1:N-1
    timeSample1     = (n-1  + epsOld )  *((periodToSample*(1 + SCO)) + timeShift);
    timeHalfSample  = (n-1/2  + eps )   *((periodToSample*(1 + SCO)) + timeShift);
    timeSample2     = (n  + eps)        *((periodToSample*(1 + SCO)) + timeShift);
   
    sample1         = interp1(windTime,windStream,timeSample1);
    halfSample      = interp1(windTime,windStream,timeHalfSample);
    sample2         = interp1(windTime,windStream,timeSample2);
    
    e(n+1)          = e(n) - 2*K*real( halfSample * (conj(sample2) - conj(sample1)) );
    eps             = e(n+1);
    epsOld          = e(n);
    outStream(n)    = sample1;
end
outStream(N)        = sample2;

% This could be moved outside of function if e is given as output parameter
% errorPlot = (e>5e-7).*e;                    % Very small errors are negligible
avg = 1/50*ones(1,50);
errorPlot = conv(e,avg,'same');
figure
plot(errorPlot)
title('Error in function of iterations (averaged)')
% To remove monotonic decrease of error in case CFO and timeShift are
% perfectly known
% diff_e = diff(e);
% figure
% plot(diff_e)
% title('diff_e')

% inStreamSampled = inStream(1:2:end);
% figure
% plot(abs(inStreamSampled-outStream))
end

