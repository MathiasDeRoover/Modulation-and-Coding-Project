function [ bitStream ] = LDPC_decoder_hard_lite( bitStream_enc, H )
%LDPC_decoder_lite Lite version of the LDPC hard decoder.
%   The hard version for the LDPC decoder.
[M_H,N_H] = size(H);
bitStream_enc_blocks    = reshape(bitStream_enc,N_H,numel(bitStream_enc)/N_H)';
bitStream_blocks        = zeros(numel(bitStream_enc)/N_H, M_H);

wait_bar = waitbar(0,'Decoding');
for i = 1:numel(bitStream_enc)/N_H
    v_nodes = bitStream_enc_blocks(i,:);
    
    iterate = true;
    while iterate
        c_nodes = mod(sum(v_nodes & H,2),2);
        temp = (sum(xor(c_nodes,v_nodes)&H)+bitStream_enc_blocks(i,:)) ./ (sum(H)+1);
        temp2 = temp == 0.5;
        v_nodes_new = and(round(temp),~temp2) + and(bitStream_enc_blocks(i,:),temp2);
        if v_nodes_new == v_nodes
            iterate = false;
        end
        v_nodes = v_nodes_new;
    end 
    bitStream_blocks(i,:) =  v_nodes(end-M_H+1:end);
    waitbar(i/(numel(bitStream_enc)/N_H),wait_bar);
end
close(wait_bar)
bitStream_blocks = bitStream_blocks';
bitStream = bitStream_blocks(:);
end

