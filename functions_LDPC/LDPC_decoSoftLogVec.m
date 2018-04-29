function [ recStream ] = LDPC_decoSoftLogVec( inStream, H, sig, it_lim )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[j_num,i_num] = size(H);
if(mod(numel(inStream),i_num)~=0)
    error('SoftLog decoding error: encoded bitstream not a multitude of block size')
end

inStream_block  = reshape(inStream,i_num,[])';

H   = reshape(H',1,i_num,[]);
c   = inStream_block;
Lq  = -2*c/(sig^2);
Lq  = repmat(Lq,1,1,j_num);

Lc = Lq;

it = 0;
while (it<it_lim) && any(any(mod(c * permute(H,[2,3,1]),2)))
    
    Lq_sign = sign(Lq);
    Lq_abs  = abs(Lq);
    
    Lr_step1 = Lq_sign;
    Lr_step1(~repmat(H,size(Lr_step1,1),1,1)) = 1;
    Lr_step2 = prod(Lr_step1,2);
    Lr_step3 = repmat(Lr_step2,1,i_num,1);
    Lr_step4 = Lr_step3./Lr_step1;
    
    Lr_step5 = log10( (exp(Lq_abs)+1) ./ (exp(Lq_abs)-1) );
    Lr_step5(~repmat(H,size(Lr_step5,1),1,1)) = 0;
    Lr_step6 = sum(Lr_step5,2);
    Lr_step7 = Lr_step6 - Lr_step5;
    Lr_step8 = log10( (exp(Lr_step7)+1) ./ (exp(Lr_step7)-1) );
    
    Lr = Lr_step4 .* Lr_step8;
    
    
    Lq_step1 = Lr;
    Lq_step1(~repmat(H,size(Lq_step1,1),1,1)) = 0;
    Lq_step2 = sum(Lq_step1,3);
    Lq_step3 = repmat(Lq_step2,1,1,j_num);
    Lq_step4 = Lq_step3 - Lq_step1;
    
    Lq = Lc + Lq_step4;
    
    
    LQ = Lc + Lq_step3;
    
    c = LQ(:,:,1)<0;
    
    it = it + 1;
end

recStream_block = c(:,end-j_num+1:end);
recStream = double(reshape(recStream_block',[],1));

end

