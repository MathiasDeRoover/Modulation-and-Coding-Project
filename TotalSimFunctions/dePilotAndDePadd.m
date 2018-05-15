function [dataOut,deltaf] = dePilotAndDePadd(dataIn,paddLength,pilotLength,dataLength,pilot,K,T,varargin)
%DEPILOTANDDEPADD Depilot and depadd
%   Depilot, depadd, frame and frequency aquisistion
dataIn = dataIn(:);
pilot = pilot(:);
T = T(:);

% ////////////////////////////
%   Find Dk(n) for all k and n
% ////////////////////////////
dataLen = numel(dataIn);
N   = numel(pilot);
K_vec   = (1:K)';
Dk  = zeros(numel(K_vec),dataLen);
for ki = 1:numel(K_vec)
    k          = K_vec(ki);
    n          = 0:dataLen-N;
    l          = (k:N-1)';
    Dk(ki,n+1) = 1/(N-k).*sum( (conj(dataIn(n+l+1)).*pilot(l+1)) .* conj(conj(dataIn(n+l-k+1)).*pilot(l-k+1)),1);
end

% ////////////////////////////
%   Find the pilots in dataIn
% ////////////////////////////
Dk_abs          = sum(abs(Dk));
Dk_block_abs    = reshape(Dk_abs,pilotLength+dataLength,[]);

[~,block_n_est] = max(Dk_block_abs);
n_est = block_n_est + (0:numel(block_n_est)-1)*(pilotLength+dataLength);

% ////////////////////////////
%   Estimate delta f using Dk
% ////////////////////////////
Dk_angle    = sum(angle(Dk)./K_vec);
deltaf      = median(-1./(2*K*pi*T(n_est)') .* Dk_angle(n_est));

% ////////////////////////////
%   Comansate the found delta f
% ////////////////////////////
disp(deltaf*2*pi)
time_vec = cumsum(T);
time_block_vec = reshape(time_vec,pilotLength+dataLength,[]);
dataIn_block = reshape(dataIn,pilotLength+dataLength,[]);
dataOut_block = dataIn_block .* exp(-1i*2*pi*deltaf*time_block_vec);
dataOut = dataOut_block(:);

% ////////////////////////////
%   linear interpolate pilot phase
% ////////////////////////////
tmp = [dataOut_block;zeros(size(dataOut_block))];
phas = zeros(1,numel(dataOut));
phas(n_est(1):n_est(1)+pilotLength-1) = angle(tmp(block_n_est(1):block_n_est(1)-1+pilotLength,1));
for i=2:numel(n_est)
    phas(n_est(i):n_est(i)+pilotLength-1) =  angle(tmp(block_n_est(i):block_n_est(i)-1+pilotLength,i))-angle(tmp(block_n_est(i-1):block_n_est(i-1)-1+pilotLength,i-1));
end
phas = interp1(n_est,phase_vec,1:numel(dataOut),'linear','extrap');
dataOut = dataOut .* exp(1i*phas');

% ////////////////////////////
%   delete extra data
% ////////////////////////////
dataOut = [dataOut;zeros(pilotLength,1)];
dataOut(n_est'+(0:numel(pilot)-1)) = [];        % Delete pilots
dataOut(end-paddLength+1-pilotLength:end) = []; % Delete padding

% ////////////////////////////
%   graph plots if specified
% ////////////////////////////
if nargin>7
    plotvar = char(varargin{1});
else
    plotvar = 'noPlot';
end
Dk_absSum = Dk_block_abs(:);
switch plotvar
    case 'plot'
        figure
        hold on
        plot(Dk_absSum)
        tmp = zeros(size(Dk_absSum));
        tmp(n_est) = Dk_absSum(n_est);
        stem(tmp)
        ylabel('Dk-mean(Dk)');
        title('Dk and estimated loc of pilots')
        hold off
end
end

