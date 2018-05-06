function [outStream] = Channel(inStream,noisePower)
%CHANNEL Sends inStream through channel resulting in outStream.
%   Sends inStream through channel resulting in outStream.
noise = sqrt(noisePower/2) * (randn(numel(inStream),1) + 1i*randn(numel(inStream),1));
outStream = inStream + noise;
end

