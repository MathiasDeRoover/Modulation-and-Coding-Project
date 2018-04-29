function [bitStream] = LDPC_decoder_log(bitStream_enc, H)
%LDPC_decoder_soft Soft decoder
%   Soft decoder for the LDPC encoding scheme.

[j_num,i_num] = size(H);
y = bitStream_enc;
sig = 0.5;

% % Calculation of p
% p = zeros(1,numel(y));
% for i = 1:i_num
%     p(i) = 1 / (1 + exp(-2 * y(i) / sig^2));
% end
% 
% Lc=log10()

% Step 1: Initialize with soft inputs from receiver
Lq = zeros(i_num,j_num);
for i = 1:i_num
    for j = 1:j_num
        Lq(i,j) = 2 * y(i) / sig / sig;
    end
end
Lc=Lq;
beta=abs(Lq);
alpha=sign(Lq);

iterate = true;
it_lim  = 200;
it      = 0;
while iterate
    it = it+1;
    
    % Step2: The check nodes calculate their response messages rji.
    % They calculate the probability that there is an even number
    % of 1’s among the variable nodes except ci.This probability
    % is equal to the probability rji(0) that ciis a 0.
    Lr = zeros(j_num,i_num);
    for i = 1:i_num
        for j = 1:j_num
            Rj      = find(H(j,:));  %Yields indices of nonzero elements
            Rji     = Rj(Rj~=i);     %Excludes current index i
            Lr(j,i) = prod(alpha(Rji,j)*min(beta(Rji,j)));
        end
    end

    
    % Step3: The variable nodes update their response messages to the
    % check nodes. Ci\j means all check nodes except fj.
    for i = 1:i_num
        for j = 1:j_num
            Ci      = find(H(:,i));
            Cij     = Ci(Ci~=j);
            %prod: product of all elements of argument
            Lq(i,j) = Lc(i,j) + sum(Lr(Cij,i)); %You could say sum(X,2)
        end
    end
    % Step4: Soft decision using pi and the rij
    for i = 1:i_num
        Ci      = find(H(:,i));
        LQ(i) = Lc(i,1) + sum(Lr(Ci,i));
    end
    % Step5: Hard decision
    c = zeros(1,i_num);
    for i = 1:i_num
        if LQ(i) <0
            c(i) = 1;
        else
            c(i) = 0;
        end
    end
    
    if( it == it_lim || ~any(c * H'))
        iterate = false;
    end
end
bitStream = c;

end


