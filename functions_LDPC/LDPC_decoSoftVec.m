function [ recStream ] = LDPC_decoSoftVec( inStream, H, sig, it_lim )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[j_num,i_num] = size(H);
if(mod(numel(inStream),i_num)~=0)
    error('Soft decoding error: encoded bitstream not a multitude of block size')
end

inStream_block  = reshape(inStream,i_num,[])';

c   = inStream_block;
H   = reshape(H',1,i_num,[]);
q0  = 1 ./ (1 + exp(2 * c ./ sig^2));
q0  = repmat(q0,1,1,j_num);
q1  = 1-q0;

p = q1;

it = 0;
while (it<it_lim) && any(any(mod(c * permute(H,[2,3,1]),2)))
    
    r0_step1 = 1-2*q1;
    r0_step1(~repmat(H,size(r0_step1,1),1,1)) = 1;
    r0_step1(r0_step1==0) = 1e-6;
    r0_step2 = prod(r0_step1,2);
    r0_step3 = repmat(r0_step2,1,i_num,1);
    r0_step4 = r0_step3./r0_step1;
    r0       = 1/2 + 1/2 * r0_step4;
    r1 = 1-r0;
    
    
    
    q1_step1 = r1;
    q1_step1(~repmat(H,size(q1_step1,1),1,1)) = 1;
    q1_step1(q1_step1==0) = 1e-6;
    q1_step2 = prod(q1_step1,3);
    q1_step3 = repmat(q1_step2,1,1,j_num);
    q1_step4 = q1_step3./q1_step1;
    
    q1_step5 = r0;
    q1_step5(~repmat(H,size(q1_step5,1),1,1)) = 1;
    q1_step5(q1_step5==0) = 1e-6;
    q1_step6 = prod(q1_step5,3);
    q1_step7 = repmat(q1_step6,1,1,j_num);
    q1_step8 = q1_step7./q1_step5;
    
    q1  = p .* q1_step4;
    q0  = (1 - p) .* q1_step8;
    q_K = (q1+q0);
    q1  = q1./q_K;
    q0  = q0./q_K;
    

    
    Q1  = p .* q1_step3;
    Q0  = (1-p) .* q1_step7;
    Q_K = (Q1+Q0);
    Q1  = Q1./Q_K;
    Q0  = Q0./Q_K;
    
    
    
    c = Q1(:,:,1)>0.5;
    
    
    it = it + 1;
end

recStream_block = c(:,end-j_num+1:end);
recStream = double(reshape(recStream_block',[],1));

end

