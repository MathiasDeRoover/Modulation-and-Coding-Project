function [outStream] = AddCFOAndPhase(inStream,streamT,CFO,phaseOffset)
%ADDCFOANDPHASE Add CFO and phase offset
%   Add CFO and phase offset
t = ((1:numel(inStream))-1)*streamT;
outStream = inStream .* exp( 1i * (CFO*t' + phaseOffset) );
end

