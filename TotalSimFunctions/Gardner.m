function [outStream] = Gardner(inStream, K)
%GARDNER Summary of this function goes here
%   Detailed explanation goes here
e(1) = 0;
for i = 1:numel(inStream)/2
    e(i+1) = e(i) + 2*K*real((inStream(i+1)+e(i))*( conj(inStream(i+2)+e(i)) - conj(inStream(i)+e(i)) ));
end

