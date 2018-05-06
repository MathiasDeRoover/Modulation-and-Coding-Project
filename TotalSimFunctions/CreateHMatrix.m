function [H0] = CreateHMatrix(c_length,v_length)
%CREATEHMATRIX Create inital parity check matrix
%   Create inital parity check matrix
H0 = makeLdpc(c_length,v_length,0,1,3);
end
