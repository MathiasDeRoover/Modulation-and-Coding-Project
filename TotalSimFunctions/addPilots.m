function [ outStream, pilotMapped ] = addPilots( symStream, pilotL, nrOfPilots, modu)
% This function will add pilots to the symbolStream which will be used to
% remove the phase offset from the data. Pilots will look like this: 
% 101010 for a pilot length of 6


% Create pilot
pilot = zeros(pilotL,1);
pilot(2:2:end) = 1;

% Transform pilot to symbols
[modulation,bps] = ModuToModulation(modu);

pilotMapped = mapping(pilot,bps,modulation);

N = numel(symStream);
dataLength = round(N/nrOfPilots);
outStream = [];
for i=0:nrOfPilots-1
    try
        data = symStream(i*dataLength+1:(i+1)*dataLength);
    catch
        data = symStream(i*dataLength+1:end);
    end
    outStream = [outStream; pilotMapped; data];    
end



% symStreamResh = reshape(symStream,nrOfPilots,[]); % Does not work if
% symStream is not divisable by nrOfPilots




% An alternative to the nrOfPilots is to use a pilot frequency.
% The pilotFreq is a value between 0 and 1 which represents how frequent
% you want pilots intertwined with data. For a value of 0.05, 5% of your
% final stream will consist of pilots.

% % Determine number of pilots based on the pilotFreq
% N = numel(symStream);
% nrOfPilotsPerSym = pilotFreq/pilotL; 
% nrOfPilots = nrOfPilotsPerSym*(N+pilotL);
% 
% % Check if we get desired pilot frequency
% totalPilotLength = nrOfPilots*pilotL;
% totalLength = totalPilotLength+N;
% totalPilotLength/totalLength            % Should be equal to pilotFreq (it is roughly equal)

end

