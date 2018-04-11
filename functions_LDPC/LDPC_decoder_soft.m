function [bitStream] = LDPC_decoder_soft(bitStream_enc, H)
%LDPC_decoder_soft Soft decoder
%   Soft decoder for the LDPC encoding scheme.

[j_num,i_num] = size(H);
y = bitStream_enc;
sig = 0.5;

% Step1: Initialize with soft inputs from receiver
q1 = zeros(i_num,j_num);
for i = 1:i_num
    for j = 1:j_num
        q1(i,j) = 1 / (1 + exp(-2 * y(i) / sig^2));
    end
end
q0 = 1-q1;
p = q1(:,1);

iterate = true;
it_lim  = 20;
it      = 0;
while iterate
    it = it+1;
    
    % Step2: The check nodes calculate their response messages rji.
    % They calculate the probability that there is an even number
    % of 1’s among the variable nodes except ci. This probability
    % is equal to the probability rji(0) that ci is a 0.
    r0 = zeros(j_num,i_num);
    for i = 1:i_num
        for j = 1:j_num
            Rj      = find(H(j,:));
            Rji     = Rj(Rj~=i);
            r0(j,i) = 1/2 + 1/2 * prod( 1 - 2 * q1(Rji,j));
        end
    end
    r1 = 1-r0;
    
    % Step3: The variable nodes update their response messages to the
    % check nodes. Ci\j means all check nodes except fj.
    K = zeros(i_num,j_num);
    for i = 1:i_num
        for j = 1:j_num
            Ci      = find(H(:,i));
            Cij     = Ci(Ci~=j);
            K(i,j)  = 1/( (1-p(i))*prod( r0(Cij,i) ) + p(i) * prod( r1(Cij,i) ) );
            q1(i,j) = K(i,j) * p(i) * prod( r1(Cij,i) );
            q0(i,j) = K(i,j) * (1-p(i)) * prod( r0(Cij,i) );
        end
    end
    
    % Step4: Soft decision using pi and the rij
    Q1 = zeros(1,i_num);
    Q0 = zeros(1,i_num);
    for i = 1:i_num
        Ci      = find(H(:,i));
        K(i)    = 1/( (1 - p(i)) * prod( r0(Ci,i) ) +  p(i) * prod( r1(Ci,i) ) );
        Q1(i)   = K(i) * p(i) * prod( r1(Ci,i) );
        Q0(i)   = K(i) * (1 - p(i)) * prod( r0(Ci,i) );
    end
    
    % Step5: Hard decision
    c = zeros(1,i_num);
    for i = 1:i_num
        if Q1(i) > 0.5
            c(i) = 1;
        else
            c(i) = 0;
        end
    end
    
    if( it == it_lim || ~any(c * H'))
        iterate = false;
    end
end

bitStream = c';

end


