function [dataOut,paddLength,pilot] = PilotAndPadd(dataIn,dataBlockLength,pilotBlockLength,modu)
%PILOTANDPADD Interweaves pilots and padds.
%   Interweaves pilots and padds

[modulation,bps] = ModuToModulation(modu);

pilot = mapping((randn(1,bps*pilotBlockLength)>0.5)',bps,modulation);
dataIn = dataIn(:);             % Convert to column

paddLength = dataBlockLength - mod(numel(dataIn),dataBlockLength);  % Determine paddlength

dataIn = [dataIn;zeros(paddLength,1)];  % padd

dataBlock   = reshape(dataIn,dataBlockLength,[]);
pilotBlock  = repmat(pilot,1,size(dataBlock,2)); 

dataOutBlock = [dataBlock;pilotBlock];
dataOut = dataOutBlock(:);
end

