function [dataOut] = dePilotAndDePadd(dataIn,paddLength,pilotLength,dataLength,pilot,K,T,varargin)
%DEPILOTANDDEPADD Depilot and depadd
%   Depilot, depadd, frame and frequency aquisistion
dataIn = dataIn(:);
pilot = pilot(:);
T = T(:);

dataLen = numel(dataIn);
N   = numel(pilot);
K_vec   = (1:K)';
Dk  = zeros(numel(K_vec),dataLen);
% for ki = 1:numel(K_vec)
%     k        = K_vec(ki);
%     n        = k+1:dataLen-N;
%     l        = (k+1:N)';
%     Dk(ki,n) = 1/(N-k).*sum((conj(dataIn(n+l-1)).*repmat(pilot(l),1,size(n,2))).*conj(conj(dataIn(n+l-1-k)).*repmat(pilot(l-k),1,size(n,2))),1);
% end
for ki = 1:numel(K_vec)
    k        = K_vec(ki);
    n        = k:dataLen-N+1;
    l        = (k:N-1)';
    Dk(ki,n) = 1/(N-k).*sum( (conj(dataIn(n+l)).*pilot(l+1)) .* conj(conj(dataIn(n+l-k)).*pilot(l-k+1)),1);
end


Dk_absSum = sum(abs(Dk));

[~,n_est] = max(Dk_absSum);
delta_f   = -1/K * sum(  angle(Dk(K_vec,n_est)) ./ (2*pi*K_vec*T(n_est)) );

time_vec  = cumsum(T);
time_vec  = [0;time_vec(1:end-1)];

dataOut = dataIn .* exp( -1i * 2*pi*delta_f*time_vec );
%dataOut = dataIn .* exp( -1i * 4000*time_vec );

% [~,n_est] = findpeaks([0,Dk_absSum],'MinPeakHeight',0.3,'MinPeakProminence',(max(Dk_absSum)-mean(Dk_absSum))*3/4);
% n_est = n_est-1;

n_est = [ fliplr(n_est:-(pilotLength+dataLength):1), (n_est+pilotLength+dataLength:pilotLength+dataLength:dataLen) ]

% 
% delta_f = -1/K .* sum( angle(Dk(K_vec,n_est)) ./ (2*pi*K_vec*T(n_est)'));
% delta_f = interp1(n_est,delta_f,1:numel(dataIn),'linear',0)';
% dataOut = dataIn .* exp( -1i*2*pi* delta_f );
%
dataOut = [dataOut;zeros(pilotLength,1)];

% angle_pilots = mean(wrapTo2Pi(angle(dataOut(n_est+(0:numel(pilot)-1)'))) - wrapTo2Pi(angle(pilot)));
% angle_data   = interp1(n_est,angle_pilots,1:numel(dataOut),'linear',0)';
% dataOut      = dataOut .* exp(-1i*angle_data);

dataOut(n_est'+(0:numel(pilot)-1)) = [];        % Delete pilots
dataOut(end-paddLength+1-pilotLength:end) = []; % Delete padding

if nargin>7
    plotvar = char(varargin{1});
else
    plotvar = 'noPlot';
end

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

