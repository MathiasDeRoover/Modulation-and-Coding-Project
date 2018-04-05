function [ bitStream ] = LDPC_decoder_hard_lite( bitStream_enc, H )
%LDPC_decoder_lite Lite version of the LDPC hard decoder.
%   The hard version for the LDPC decoder.
[M_H,N_H] = size(H);
bitStream_enc_blocks    = reshape(bitStream_enc,N_H,numel(bitStream_enc)/N_H)';
bitStream_blocks        = zeros(numel(bitStream_enc)/N_H, M_H);

for i = 1:numel(bitStream_enc)/N_H
    v_nodes = bitStream_enc_blocks(i,:);
    
    iterate = true;
    while iterate
        c_nodes = mod(sum(v_nodes & H,2),2);
        v_nodes_new = round(((sum(and(xor(c_nodes,v_nodes),H)) + bitStream_enc_blocks(i,:)) ./ (sum(H)+1))-0.001);
        %v_nodes_new = round((sum(and(xor(c_nodes,v_nodes),H))) ./ (sum(H)));
        if v_nodes_new == v_nodes
            iterate = false;
        end
        v_nodes = v_nodes_new;
    end 
    bitStream_blocks(i,:) =  v_nodes(end-M_H+1:end);
end
bitStream_blocks = bitStream_blocks';
bitStream = bitStream_blocks(:);
end

