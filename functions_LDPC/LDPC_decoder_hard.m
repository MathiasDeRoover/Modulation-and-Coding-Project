function [ bitStream ] = LDPC_decoder_hard( bitStream_enc, H )
%LDPC_DECODER_HARD Decode LDPC encoding.
%   Decode LDPC encoding using a hard decoding scheme. This method should
%   go fast and be pretty optimized.
[c_num,v_num]       = size(H);
bitstrm_enc_rshp    = reshape(bitStream_enc,v_num,[])';
H                   = reshape(H',1,v_num,[]);

v_nodes             = bitstrm_enc_rshp;
v_nodes_old         = inf(size(v_nodes));
while_it            = 0;
while_it_limit      = 10;
while (while_it ~= while_it_limit) && any(any(v_nodes - v_nodes_old))
    while_it        = while_it + 1;
    c_nodes         = mod(sum(v_nodes & H,2),2);
    c_nodes         = reshape(c_nodes,[],1,c_num);
    
    c_xor_v         = xor( v_nodes , c_nodes );
    c_xor_v_masked  = c_xor_v & H;
    avarage         = ( sum(c_xor_v_masked,3) + bitstrm_enc_rshp ) ./ ( sum(H,3) + 1 );
    avarage_mask    = avarage==0.5;
    v_nodes_old     = v_nodes;
    v_nodes         = and(round(avarage),~avarage_mask) + and(bitstrm_enc_rshp,avarage_mask);
end
bitStream           = reshape(v_nodes(:,end-c_num+1:end)',1,[])';
end

