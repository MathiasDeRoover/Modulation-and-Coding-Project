function [ bitStream_enc , newH] = LDPC_encoder_lite( bitStream, H )
%LDPC_encoder_lite LDPC encoder for small size blocks
%   LDPC encoder for small size blocks. Useful to study the main properties
%   of LDPC encoders.

%% Initializing
H         = mod(rref(H),2);         % Write H as H = [I P'] (using reduced row echelon form)
[M,N]     = size(H);
P         = H(:,N-M+1:end);         % Find P
P         = P';                     %
[M_P,~]   = size(P);

G         = [P eye(M_P)];           % Create the generator matrix

%% Encoding the bitStream
if mod(numel(bitStream),M_P) ~= 0   % Make sure that bitStream has a multitude of bits compared to M_P
    error(['Number of bits not a multiple of ',num2str(M_P)]);
end

L = numel(bitStream)/M_P;
bitStream_block = reshape(bitStream',M_P,L)';
bitStream_enc_block = zeros(L,N);
for i = 1:L
    bitStream_enc_block(i,:) = mod(bitStream_block(i,:)*G,2);
end
bitStream_enc_block = bitStream_enc_block';
bitStream_enc = bitStream_enc_block(:); 
newH = H;
end