clear
close all
clc
addpath(genpath(pwd))
%% INITIALIZATION %%
modu = 'QPSK';
[modulation,bps]    = ModuToModulation(modu);
ftaps               = 200;              % Amount of causal (and non causal) filter taps
symRate             = 2*1e6;            % 2 * cutoffFrequency = 1/T (We use the -3dB point as cutoffFrequency)
T                   = 1/symRate;        % Symbol period
M                   = 100;              % UpSample factor
fs                  = symRate * M;      % Sample Frequency
beta                = 0.3;              % Roll off factor
N                   = 128*12;           % Amount of bits in original stream
hardDecodeIter      = 10;               % Iteration limit for hard decoder
cLength             = 128;
vLength             = 256;

fc                  = 2e9;              % 2 GHz carrier freq is given as example in the slides
ppm                 = 10e-6;            % Parts per million
deltaW10            = fc*ppm;           % Carrier frequency offset CFO 10ppm 10e-6
deltaW50            = 5*deltaW10;
deltaW0             = 0;
phi0                = 0;                % Phase offset
delta               = 0;                % Sample clock offset SCO
t0                  = 0.5*T;            % Time shift
K_gard              = 0.01;             % K for Gardner

SNRdB               = 200;              % Signal to noise in dB
nrOfIterations      = 20;               % To take mean and std


for ii=1:nrOfIterations
    %% Create bitstream
    [stream_bit]        = CreateBitStream(N,1);

    %% Mapping
    stream_mapped       = mapping(stream_bit, bps, modulation);

    %% Upsample
    stream_upSampled    = upsample(stream_mapped,M);

    %% Create window
    [g,g_min]           = CreateWindow(T, fs, ftaps, beta);

    %% Apply window
    stream_wind         = conv(stream_upSampled,g);

    %% Sending through channel
    noisePower          = CalcNoisePower(stream_bit, stream_wind,1/fs,SNRdB,fs,ftaps);
    stream_channel      = Channel(stream_wind,noisePower);
%         stream_channel = stream_wind; %%%%!!!!!%%%%

    %% CFO + Phase offset
    stream_rec_CFOPhase0 = AddCFOAndPhase(stream_channel,fs,deltaW0,phi0);
    stream_rec_CFOPhase10 = AddCFOAndPhase(stream_channel,fs,deltaW10,phi0);
    stream_rec_CFOPhase50 = AddCFOAndPhase(stream_channel,fs,deltaW50,phi0);

    %% Apply window
    stream_rec_wind0     = conv(stream_rec_CFOPhase0,g_min);
    stream_rec_wind10     = conv(stream_rec_CFOPhase10,g_min);
    stream_rec_wind50     = conv(stream_rec_CFOPhase50,g_min);

    %% Truncate & Sample + Sample clock offset and time shift
    stream_rec_sample   = TruncateAndSample(stream_rec_wind0,ftaps,T,fs,delta,t0);

    %% Gardner
    disp(['Gardner iteration: ',num2str(ii)])
    [~, ~, epsilon(ii,:)]  = Gardner(stream_rec_sample, K_gard, stream_rec_wind0, ftaps, fs, T,delta, t0);
    [~, ~, epsilon2(ii,:)]  = Gardner(stream_rec_sample, K_gard, stream_rec_wind10, ftaps, fs, T,delta, t0);
    [~, ~, epsilon3(ii,:)]  = Gardner(stream_rec_sample, K_gard, stream_rec_wind50, ftaps, fs, T,delta, t0);
end

error = (t0/T-abs(epsilon))*T;
error2 = (t0/T-abs(epsilon2))*T;
error3 = (t0/T-abs(epsilon3))*T;


avg = mean(error,1);
stddev = std(error,1);
avg2 = mean(error2,1);
stddev2 = std(error2,1);
avg3 = mean(error3,1);
stddev3 = std(error3,1);


steps = 1:20:numel(avg);
figure
plot(steps,avg(steps),'-or')
hold on
plot(steps,avg2(steps),'-xb')
plot(steps,avg3(steps),'-^k')

plot(steps,avg(steps)-stddev(steps),'--or')
plot(steps,avg(steps)+stddev(steps),'--or')
plot(steps,avg2(steps)-stddev2(steps),'--xb')
plot(steps,avg2(steps)+stddev2(steps),'--xb')
plot(steps,avg3(steps)-stddev3(steps),'--^k')
plot(steps,avg3(steps)+stddev3(steps),'--^k')

title('QPSK, 2 Mbps symbol rate, 0.3 roll-off, no noise')
xlabel('Symbols')
ylabel('Time error (mean \pm stdv)')
ylim(1e-7*[-0.5,3])
legend('no CFO','CFO = 10 ppm','CFO = 50 ppm')
hold off

rmpath(genpath(pwd))