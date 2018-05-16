function [dataOut,deltaf,n_est] = removePhase(dataIn,paddLength,pilotLength,dataLength,pilot,K,T,varargin)
%DEPILOTANDDEPADD Depilot and depadd
%   Depilot, depadd, frame and frequency acquisistion
dataIn = dataIn(:);
pilot = pilot(:);
T = T(:);

% ////////////////////////////
%   Find Dk(n) for all k and n
% ////////////////////////////
dataLen = numel(dataIn);
N   = numel(pilot);
K_vec   = (1:K).';
Dk  = zeros(numel(K_vec),dataLen);
for ki = 1:numel(K_vec)
    k          = K_vec(ki);
    n          = 0:dataLen-N;
    l          = (k:N-1).';
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
deltaf      = median(-1./(2*K*pi*T(n_est).') .* Dk_angle(n_est));

% ////////////////////////////
%   Comansate the found delta f
% ////////////////////////////
time_vec = cumsum(T);
time_block_vec = reshape(time_vec,pilotLength+dataLength,[]);
dataIn_block = reshape(dataIn,pilotLength+dataLength,[]);
dataOut_block = dataIn_block .* exp(-1i*2*pi*deltaf*time_block_vec);
dataOut_block_copy = dataIn_block;
dataOut = dataOut_block(:);
dataOut_copy = dataOut_block_copy(:);


% ////////////////////////////
%   linear interpolate pilot phase
% ////////////////////////////
% tmp = [dataOut_block;zeros(size(dataOut_block))];
% phase_vec = zeros(1,numel(n_est));
% phase_vec(1) = 0;
% for i=2:numel(n_est)
%     phase_vec(i) = mean(angle(tmp(block_n_est(i):block_n_est(i)-1+pilotLength,i))-angle(tmp(block_n_est(i-1):block_n_est(i-1)-1+pilotLength,i-1)));
% end
% phas = interp1(n_est,phase_vec,1:numel(dataOut),.'linear.',.'extrap.');

dataOut = [dataOut;zeros(pilotLength,1)];
dataOut_copy = [dataOut_copy;zeros(pilotLength,1)];


tmp_n = n_est.'+(0:numel(pilot)-1);
tmp = [pilot.';dataOut(tmp_n)];
tmp1 = unwrap((angle(tmp)));
tmp_diff = diff(tmp1);
%tmp2 = (angle(pilot.'));
%tmp_ang = angle( exp(1i*tmp1).*exp(-1i*tmp2));
tmp_ang = cumsum(tmp_diff);
rest_angle = tmp_ang(:);

phas = interp1(tmp_n(:),rest_angle(:),1:numel(dataOut),'linear',0);

dataOutTmp = dataOut;
% dataOut = dataOut .* exp(-1i*phas.');
dataOut = dataOut_copy .* exp(-1i*phas.');
dataOut = dataOut_copy;


% ////////////////////////////
%   delete extra data
% ////////////////////////////

scatter_1 = [pilot.';dataOut(n_est.'+(0:numel(pilot)-1))];
scatter_2 = [pilot.';dataOutTmp(n_est.'+(0:numel(pilot)-1))];

dataOut(n_est.'+(0:numel(pilot)-1)) = [];        % Delete pilots
dataOut(end-paddLength+1-pilotLength:end) = []; % Delete padding
end

