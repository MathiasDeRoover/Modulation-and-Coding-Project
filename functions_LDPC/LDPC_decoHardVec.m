function [ bitStream ] = LDPC_decoHardVec( bitStream_enc, H ,while_it_limit )
%LDPC_DECODER_HARD Decode LDPC encoding.
%   Decode LDPC encoding using a hard decoding scheme. This method should
%   go fast and be pretty optimized.
[c_num,v_num]       = size(H);
bitstrm_enc_rshp    = reshape(bitStream_enc,v_num,[])';
H                   = reshape(H',1,v_num,[]);     

v_nodes             = bitstrm_enc_rshp;
v_nodes_change      = ones(size(v_nodes))*0.5;
while_it            = 0;

if ~exist('while_it_limit','var')
        while_it_limit      = 10;
end

while (while_it < while_it_limit) && any(any(mod(v_nodes * permute(H,[2,3,1]),2)))
    while_it        = while_it + 1;
    c_nodes         = mod(sum(v_nodes & H,2),2);
    c_nodes         = reshape(c_nodes,[],1,c_num);
    
    c_xor_v         = xor( v_nodes , c_nodes );
    c_xor_v_masked  = c_xor_v & H;
    average         = ( sum(c_xor_v_masked,3) + v_nodes + v_nodes_change ) ./ ( sum(H,3) + 2 );
    average_mask    = average == 0.5;
    v_nodes_old     = v_nodes;
    v_nodes         = and(round(average),~average_mask) + and(v_nodes,average_mask);
    v_nodes_change = v_nodes;
    v_nodes_change(v_nodes_old == v_nodes) = ~v_nodes(v_nodes_old == v_nodes);
    
end
bitStream           = reshape(v_nodes(:,end-c_num+1:end)',1,[])';
end

