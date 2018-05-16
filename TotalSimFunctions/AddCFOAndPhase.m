function [outStream] = AddCFOAndPhase(inStream,streamFrequency,CFO,phaseOffset)
%ADDCFOANDPHASE Add CFO and phase offset
%   Add CFO and phase offset
t = ((1:numel(inStream))-1)/streamFrequency;
outStream = inStream .* exp( 1i * (2*pi*CFO*t' + phaseOffset) );
end

