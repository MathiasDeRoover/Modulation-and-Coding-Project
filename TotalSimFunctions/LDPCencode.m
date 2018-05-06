function [outStream,newH] = LDPCencode(inStream,H0)
%LDPCENCODE Encode inStream with H0
%   Encode inStream with H0
bitStream_blk               = reshape(inStream,size(H0,1),[]);
[bitStream_cod_blk,newH]    = makeParityChk(bitStream_blk, H0, 0);  % Create parity check bits and reshape H
bitStream_cod_blk           = [bitStream_cod_blk;bitStream_blk];    % Unite parity check bits and message
outStream                   = reshape(bitStream_cod_blk,[],1);
end

