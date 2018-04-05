function [bitStream] = CreateBitStream(N,bps)
%UNTITLED4 Creates a random bitStream of N symbols with bps symbols per bit
%   Creates a random bitStream of N symbols with bps symbols per bit
bitStream = double(rand(N*bps,1) > 0.5);
end

