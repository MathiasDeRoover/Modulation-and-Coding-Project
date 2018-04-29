function [ received ] = hardDecoderLDPC_G2( received, H, messageBlockSize, codedBlockSize )
%hardDecoderLPDC Summary of this function goes here
%   Receives a 1 x (nmbrMessages.codedBlockSize) matrix and gives also a
%   matrix of this format back. For the calculations, it is reshaped
%   internally to a nmbrMessages x codedBlockSize sized matrix

received = reshape(received,codedBlockSize,[])';

nmbrMessages = size(received,1);
maxIterations = 150;
actualIteration = 0;

decisionvector = ones(nmbrMessages,codedBlockSize);

while actualIteration<maxIterations
    newreceived = received; % v_nodes; variable nodes; received codeword; c's
    %decisionvector = ones(nmbrMessages,codedBlockSize); %If positive, keep the bit of that element, if negative, change it, if zero (no majority), also keep
    for i = 1:messageBlockSize % i goes through all columns of H -> through all check nodes
        parity = mod(newreceived * H(i,:)',2); % XOR of c's corresponding to ones in H - c_nodes; check nodes; fi
        Hfull = repmat(H(i,:),[nmbrMessages,1]);
        
        if size(find(parity),1)~=0
            decisionvector(find(parity),:) = decisionvector(find(parity),:) - Hfull(find(parity),:);
        end
        if size(find(~parity),1)~=0
            decisionvector(find(~parity),:) = decisionvector(find(~parity),:) + Hfull(find(~parity),:);
        end
    end
    
    decisionvector = sign(decisionvector);
    zeroindices = find(decisionvector==0);
    decisionvector(zeroindices) = ones(length(zeroindices),1);
    negindices = find(decisionvector==-1);
    decisionvector(negindices) = zeros(length(negindices),1); 
    decisionvector = not(decisionvector); % 0 = keep; 1 = change
    newreceived = xor(newreceived,decisionvector);
    
    decisionvector = decisionvector*2; % intial conditions next iteration; back to voting
    
    if newreceived == received
        break;
    else
        received = newreceived;
        actualIteration=actualIteration+1
        if actualIteration==maxIterations
            warning('hardDecoderLPDC: Max looping exceeded')
            break;
        end
    end
    
end

received = received(:,codedBlockSize-messageBlockSize+1:end);
received = reshape(received',[],1);

end

