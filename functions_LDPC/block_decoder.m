function [ bitStream ] = block_decoder( bitStream_enc, H  )
%block_decoder Block Decoder
%   Non iteraitve block decoder

[c,v] = size(H);
if(mod(numel(bitStream_enc),v)~=0)
    error('Block Decoding error: Bitstream not multitude of block size')
end

r = reshape(bitStream_enc(:),v,[])';
s = mod(r * H',2);

e = (H\s')';






r = mod(r + e,2);
bitStream_blocks = r(end-c+1:end,:)';
bitStream = bitStream_blocks(:);
end

