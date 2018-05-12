function [dataOut,paddLength,pilot] = PilotAndPadd(dataIn,dataBlockLength,pilotSymbols,pilotBlockLength,varargin)
%PILOTANDPADD Interweaves pilots and padds.
%   Interweaves pilots and padds
if(mod(pilotBlockLength,numel(pilotSymbols))~=0)
    error('pilotBlockLength needs to be a multitude of pilotSymbols length')
end

if nargin>4
   modu = varargin{1};
   [modulation,bps] = ModuToModulation(modu);
   if mod(numel(pilotSymbols),bps) ~= 0
      error(['pilotsymbol length needs to be a multitide of bps:',num2str(bps)]);
   end
   pilotSymbols = mapping(pilotSymbols',bps,modulation);
end

dataIn = dataIn(:);             % Convert to column
pilotSymbols = pilotSymbols(:); % Convert to column

paddLength = dataBlockLength - mod(numel(dataIn),dataBlockLength);  % Determine paddlength

dataIn = [dataIn;zeros(paddLength,1)];  % padd

dataBlock   = reshape(dataIn,dataBlockLength,[]);
pilotBlock  = repmat(pilotSymbols,pilotBlockLength/numel(pilotSymbols),size(dataBlock,2)); 
pilot       = pilotBlock(:,1);

dataOutBlock = [dataBlock;pilotBlock];
dataOut = dataOutBlock(:);
end

