function [bitStream] = LDPC_decoder_soft_log_BPSK(bitStream_enc, H,sig)
%LDPC_decoder_soft_log_BPSK Soft decoder in log domain
%   Soft decoder for the LDPC encoding scheme in log domain. Specifically for BPSK.
[j_num,i_num] = size(H);
if(mod(numel(bitStream_enc),i_num)~=0)
    error('Soft decoding error: encoded bitstream not a multitude of block size')
end

bitStream = zeros(numel(bitStream_enc)/i_num*j_num,1);
for bstr = 1:numel(bitStream_enc)/i_num
    y = bitStream_enc( (1:i_num) + (bstr-1)*i_num );

    %% Initialize
    L_c = - 2*y/(sig^2);
    L_q = repmat(L_c,1,j_num);
    
    iterate = true;
    it_lim  = 10;
    it      = 0;
    while iterate
        
        chi = sign(L_q);
        alpha = abs(L_q);
        
        L_r = zeros(j_num,i_num);
        for i = 1:i_num
            for j = 1:j_num
                Rj      = find(H(j,:)); %% Maybe change this
                Rji     = Rj(Rj~=i);
                %L_r(j,i) = prod( chi(Rji,j) ) * -log10(tanh( 1/2 * sum( -log10(tanh(1/2 * alpha(Rji,j)))))); %% Possible change to log10
                L_r(j,i) = prod( chi(Rji,j) ) * min( alpha(Rji,j) ); %% Possible change to log10
            end
        end
        
        L_Q = zeros(i_num,1);
        for i = 1:i_num
            for j = 1:j_num
                Ci      = find(H(:,i)); %% Maybe change this
                Cij     = Ci(Ci~=j);
                L_q(i,j) = L_c(i) + sum( L_r(Cij,i) );
            end
            L_Q(i) = L_c(i) + sum( L_r(Ci,i) );
        end
        
        c = (L_Q < 0)'; % 0= log10(1) so >1--> >0,so <1--> <0

        if( it == it_lim || ~any(mod(c * H',2)))
            iterate = false;
        end
        it = it+1;
    end
    
    bitStream((1:j_num) + (bstr-1)*j_num ) = c(end-j_num+1:end)'; %c(1:end)' ?
    % (1:j_num) + (bstr-1)*j_num :: Create right indices
    % end-j_num+1:end :: save last determined values
end
end

