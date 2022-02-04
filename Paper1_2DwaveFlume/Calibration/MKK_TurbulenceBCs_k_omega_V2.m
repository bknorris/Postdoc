%% Calculate Turbulence BCs for IHFOAM irregularWaveModel
% The model requires an estimate of TKE (k), specific dissipation rate 
% (eps), and the turbulent eddy viscosity (nut)
%
% NOTE: We don't have velocity measurements near the inlet of the model
% domain, so these estimates are derived from the Aquadopps within the HR
% and LR patches. 
%
% Inputs: waveStats files from HR and LR Aquadopps
% Outputs: turbulent boundary conditions for the model
%
% Updates: 01/28/21: Added functionality to average data over multiple
% bursts
%
% This is version 2 of this script.
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

%% Copy Bursts from MKK_Inputs_for_irregularWaveModel.m
RBRbursts = [151 158]; %RBR bursts to analyze (range: 141 - 162)
AQDbursts = [61 68]; %ADCP bursts to analyze (range: 51 - 72)

%% Calculate BCs
%burstNum is determined by C15, which is different than the AQDPs
%Use burstNum to figure out which AQDP burst to use
idxHR = find(wavesHR.wave.burst >= AQDbursts(1) & wavesHR.wave.burst <= AQDbursts(2));
idxLR = find(wavesLR.wave.burst >= AQDbursts(1) & wavesLR.wave.burst <= AQDbursts(2));

%Extract epsilon values from waveStats. Note, discard lower 3 bins (these
%are obscured by the bed)
epsilonAvg = mean((mean(wavesHR.wave.epsilon(1:11,idxHR),2)+mean(wavesLR.wave.epsilon(1:11,idxLR),2))./2);

%Compute k from epsilon (see Pope et al., 2006)
zuv = median((wavesLR.wave.zuv+wavesHR.wave.zuv)./2);
c0 = 0.19; 
kappa = 0.41;
k = ((epsilonAvg^(2/3))*(kappa*zuv^(2/3)))/c0;

%compute turbulenceLengthScale (l)
l = (k^(3/2))/epsilonAvg;

%Compute nut
Cmu=0.09;
omega = sqrt(k)/((Cmu^(1/4))*l);

nut = k/omega;

X = ['k= ',num2str(k)];
display(X);

X = ['omega= ',num2str(omega)];
display(X);

X = ['nut= ',num2str(nut)];
display(X);