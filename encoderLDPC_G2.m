function [ encodedpackets, H ] = encoderLDPC( bitsTx, messageBlockSize, codedBlockSize )
%encoderLPDC Summary of this function goes here
%   Receives a 1 x (nmbrMessages.codedBlockSize) matrix and gives also a
%   matrix of this format back

if mod(bitsTx,messageBlockSize)~=0
    error('The length of bitsTx is not a multiple of blockSize')
end

%Check if int if

H0 = makeLdpc(messageBlockSize,codedBlockSize,0,1,3); % Create initial parity check matrix of size 128 x 256

infobits = reshape(bitsTx,messageBlockSize,[]);

[codedbits, H] = makeParityChk(infobits, H0, 0); % Encode information bits (128 x number of packets) and generate final parity check matrix

encodedpackets = reshape([codedbits;infobits],[],1);

end

