function [ bitStream ] = LDPC_decoHard( bitStream_enc, H, while_it_limit )
%LDPC_decoder_lite Lite version of the LDPC hard decoder.
%   The hard version for the LDPC decoder.
[M_H,N_H] = size(H);
bitStream_enc_blocks    = reshape(bitStream_enc,N_H,numel(bitStream_enc)/N_H)';             % Reshape bitstream back to blocks of code
bitStream_blocks        = zeros(numel(bitStream_enc)/N_H, M_H);                             % Initialize vector for decoded code blocks

if ~exist('while_it_limit','var')
        while_it_limit      = 10;
end

for i = 1:numel(bitStream_enc)/N_H
    v_nodes         = bitStream_enc_blocks(i,:);                                            % Decode one block of bits at a time
    v_nodes_change  = ones(size(v_nodes)) * 0.5;
    iterations = 0;
    while iterations < while_it_limit && any(mod(v_nodes * H',2))
        c_nodes = mod(sum(v_nodes & H,2),2);                                                % Check comment below
        temp = (sum(xor(c_nodes,v_nodes)&H) + v_nodes + v_nodes_change) ./ (sum(H)+2);       % Check comment below
        temp2 = temp == 0.5;
        v_nodes_new = and(round(temp),~temp2) + and(v_nodes,temp2);
        v_nodes_change = v_nodes;
        v_nodes_change(v_nodes_old == v_nodes) = ~v_nodes(v_nodes_old == v_nodes);
        v_nodes = v_nodes_new;
        iterations = iterations + 1;
    end 
    bitStream_blocks(i,:) =  v_nodes(end-M_H+1:end);
end
bitStream_blocks = bitStream_blocks';
bitStream = bitStream_blocks(:);
end

%% Comments
% c_nodes = mod(sum(v_nodes & H,2),2);
% H contains the connections between the c_nodes and the v_nodes. Perform a
% logical & operation with H and the v_nodes to know which bits from the
% v_nodes contribute to the c_nodes.
% If you now sum over the rows of H, you will count the number of ones
% received by the c_nodes. Since you want to perform an exclusive or
% operation between all bits received at the c_nodes, an odd number of ones
% will result in a one and an even number of bits result in a zero; hence
% the modulo 2 operation
%
% temp = (sum(xor(c_nodes,v_nodes)&H)+bitStream_enc_blocks(i,:)) ./ (sum(H)+1);
% Now you want to send back the result from the check nodes to the variable
% nodes, but you want to exclude the influence of the v_node you are
% sending to (xor at c_nodes should be performed with bits from all other
% v_nodes, except the one you are sending back to). To remove this
% influence you perform an xor with the c_nodes and v_nodes.
% Parity check matrix H still contains the connections between c_nodes and
% v_nodes. By again performing a logical & between our result and H, we can
% find which c_nodes send a one back to the v_nodes. This time, summing
% over the columns gives the amount of ones received by the v_nodes.
% Next we add the initially received bits from the channel. 
% Assume there are no errors. If v1 receives a 1 from the channel and is
% connected to 3 check nodes, then it should receive 3 ones back from the
% check nodes? The received sum should thus be 3. If we add the one that
% was initially received the sum becomes 4. Dividing this by sum(H) + 1
% will again give 1 (sum(H) gives amount of connected nodes; +1 is for
% counting bit received from channel).
% IMPORTANT: doesn't this means we should use the original H instead of
% newH?
% If there is noise, the result will be something between [1,0]. By
% rounding the result we apply the majority rule and select the most
% likely result. In case the result is 0.5 we believe the bit received from
% the channel is the correct one.


