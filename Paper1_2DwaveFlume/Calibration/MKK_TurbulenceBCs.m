%% Calculate Turbulence BCs for IHFOAM irregularWaveModel
% The model requires an estimate of TKE (k), dissipation rate (eps), and
% the turbulent eddy viscosity (nut)
%
% NOTE: We don't have velocity measurements near the inlet of the model
% domain, so these estimates are derived from the Aquadopps within the HR
% and LR patches. 
%
% Inputs: waveStats files from HR and LR Aquadopps
% Outputs: turbulent boundary conditions for the model
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all;
dataPath = 'C:\Users\user\Documents\Data\InstrumentData\WavesData\';
wavesHR = load([dataPath 'MKK18HR101aqdHR_waves.mat']);
wavesLR = load([dataPath 'MKK18LR101aqdHR_waves.mat']);
C15 = load([dataPath 'MKK18C1501rbr_waves.mat']);

%burstNum is determined by C15, which is different than the AQDPs
%Use burstNum to figure out which AQDP burst to use
burstNum = 160;
burstID = find(C15.wave.burst == burstNum);
[~,timeIDHR] = min(abs(wavesHR.wave.time-C15.wave.time(burstID)));
[~,timeIDLR] = min(abs(wavesLR.wave.time-C15.wave.time(burstID)));

%Extract epsilon values from waveStats. Note, discard lower 3 bins (these
%are obscured by the bed)
epsilonAvg = mean((wavesHR.wave.epsilon(1:11,timeIDHR)+wavesLR.wave.epsilon(1:11,timeIDLR))./2);

%Compute k from epsilon (see Pope et al., 2006)
zuv = median((wavesLR.wave.zuv+wavesHR.wave.zuv)./2);
c0 = 0.19; 
kappa = 0.41;
k = ((epsilonAvg^(2/3))*(kappa*zuv^(2/3)))/c0;

%Compute nut
Cmu=0.09;
nut = Cmu * k^2 / epsilonAvg;

X = ['k= ',num2str(k)];
display(X);

X = ['epsilon= ',num2str(epsilonAvg)];
display(X);

X = ['nut= ',num2str(nut)];
display(X);