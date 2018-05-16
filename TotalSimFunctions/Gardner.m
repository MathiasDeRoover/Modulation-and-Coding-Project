function [outStream, T_est, e] = Gardner(inStream, K, windStream, ftaps, windFrequency, periodToSample,SCO, timeShift)
%GARDNER Summary of this function goes here
%   Detailed explanation goes here
% sample 161 of windStream is the first symbol that is tranmitted  
N = numel(inStream);
e = zeros(N,1);
windStream  = windStream(2*ftaps+1:end); %Comment for test

windTime    = ((1:numel(windStream))-1)/windFrequency;
epsOld  = 0;
eps     = 0;
outStream   = zeros(N,1);
T_est       = zeros(N,1);

% figure
% stem(windTime,real(windStream))
% hold on
% % plot lines at -1 and 1
% plot(windTime,ones(size(windTime)),'--')
% plot(windTime,-ones(size(windTime)),'--')
% hold off



for n = 1:N-1
    timeSample1     = (n-1  - epsOld )*(periodToSample*(1+SCO)) + timeShift;
    timeHalfSample  = (n-1/2  - eps )*(periodToSample*(1+SCO)) + timeShift;
    timeSample2     = (n  - eps )*(periodToSample*(1+SCO)) + timeShift;
   
    sampleTimes             = [timeSample1,timeHalfSample,timeSample2];
    samples                 = interp1(windTime,windStream,sampleTimes);
    sample1                 = samples(1);
    halfSample              = samples(2);
    sample2                 = samples(3);
%     sample1         = interp1(windTime,windStream,timeSample1);
%     halfSample      = interp1(windTime,windStream,timeHalfSample);
%     sample2         = interp1(windTime,windStream,timeSample2);
    
    sample1(isnan(sample1)) = 0;
    halfSample(isnan(halfSample)) = 0;
    sample2(isnan(sample2)) = 0;
    
    
    e(n+1)          = e(n) + 2*K*real( halfSample * (conj(sample2) - conj(sample1)) );
    eps             = e(n+1);
    epsOld          = e(n);
    outStream(n)    = sample1;
    T_est(n)        = timeSample2-timeSample1;
end
T_est(N)            = timeSample2-timeSample1;
outStream(N)        = sample2;

end

