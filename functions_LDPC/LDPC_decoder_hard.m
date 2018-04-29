function [ bitStream ] = LDPC_decoder_hard( bitStream_enc, H ,while_it_limit )
%LDPC_DECODER_HARD Decode LDPC encoding.
%   Decode LDPC encoding using a hard decoding scheme. This method should
%   go fast and be pretty optimized.
[c_num,v_num]       = size(H);

%%% This code is only needed if receiver does not know amount of
%%% zerosToPad. I am not 100% sure that we don't remove useful information
%%% when using this code.

% %%%
% zerosToPad = c_num - mod(length(bitStream_enc),c_num);
% if (zerosToPad ~= 0)
% %     disp('We need to pad some zeros')
%     bitStream_enc = [bitStream_enc; zeros(zerosToPad,1)];
%     while mod(numel(bitStream_enc),bps) ~= 0                % Infinite loops should not happen, but be careful!
%         bitStream_enc = [bitStream_enc; zeros(c_num,1)];
%         zerosToPad = zerosToPad + c_num;
%     end
% end
%%%

bitstrm_enc_rshp    = reshape(bitStream_enc,v_num,[])';
H                   = reshape(H',1,v_num,[]);                               % Why do you make dimensions 1x256x128?

v_nodes             = bitstrm_enc_rshp;
v_nodes_old         = inf(size(v_nodes));
while_it            = 0;

if ~exist('while_it_limit','var')
        while_it_limit      = 10;
end

while (while_it < while_it_limit) && any(any(mod(v_nodes * permute(H,[2,3,1]),2)))
    while_it        = while_it + 1;
    c_nodes         = mod(sum(v_nodes & H,2),2); 
    %v is compared to H in a logical AND way. Thereafter for each row the
    ...sum is taken and that result is modulo 2 divided to get the modulo-2 sum of each row
    c_nodes         = reshape(c_nodes,[],1,c_num); %Reshape C to be 3D array of columns
    
    %More comments: see decoder_hard_lite
    c_xor_v         = xor( v_nodes , c_nodes ); 
    c_xor_v_masked  = c_xor_v & H; 
    % Average calculates the value that would be given to the node. If 1 or
    % zero--> correct, else error. 
    average         = ( sum(c_xor_v_masked,3) + bitstrm_enc_rshp ) ./ ( sum(H,3) + 1 );
    average_mask    = average==0.5;% One where avg = 0.5
    % If average = 0.5 you assume initial value was correct (so put in old
    % value (second and below) and neglect averaged value( first and with
    % complement of mask).
    v_nodes_old     = v_nodes;
    v_nodes         = and(round(average),~average_mask) + and(bitstrm_enc_rshp,average_mask);
    % Fist term: if average ~= 0.5, use rounded average (>0.5 -->1 , less
    % than 0.5--> 0, this is the majority voting).
    % Second term: if average = 0.5 than first and is zero and second and
    % gives original value of that v_node.
end
bitStream           = reshape(v_nodes(:,end-c_num+1:end)',1,[])';
% bitStream           = bitStream(1:end-zerosToPad);
end

