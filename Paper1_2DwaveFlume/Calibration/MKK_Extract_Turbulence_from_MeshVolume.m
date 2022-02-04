%Extract modeled turbulence values from the 'Turbulence Refinement Models'.
%This script iteratively loads the data from each models (user specified)
%into a data structure and optionally plots the results as a time series
%(image). Results can also be optionally saved to disk. 
%
% Updates:
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
modelName = {'MKK_CombinedHR_level1'};
modelFigName = 'Level1 refinement model [0.3 - 0.15 m]';
% modelFigName = 'Level2 refinement model [0.3 - 0.075 m]';
% modelFigName = 'Level3 refinement model [0.3 - 0.0375 m]';
% modelFigName = 'Level4 refinement model [0.3 - 0.0187 m]';
% modelFigName = 'Level7 refinement model [0.3 - 0.002 m]';
RBRbursts = [151 158]; %RBR bursts to analyze (range: 141 - 162)
rampTime = 20; %timeseries will be adjusted by this number of seconds to skip data during 'rampTime'
modelTimeStep = 0.125; %model timestep in seconds
modelDuration = 120; %length of timeseries in seconds
saveFigs = 1; %save figures to disk; 0 is off, 1 is on
saveData = 1; %save out profile data to disk; 0 is off, 1 is on

%% Load and process the model data
%Get data paths
if length(modelName) ~= 1
    error('Only process one model at a time!')
end
if length(RBRbursts) == 1
    burstTxt = sprintf('%d',RBRbursts(1));
else
    burstTxt = sprintf('%d-%d',RBRbursts(1),RBRbursts(2));
end
modelBasePath = ['c:\Users\user\Documents\Models\MKK_Combined\2D\Burst' burstTxt '\TurbulenceRefinementModels\' modelName{:} '\MKKmodel\postProcessing\bathySample\surface\'];
modelFolders = dir(modelBasePath);
modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
modelFolders = {modelFolders.name};
sortModels = str2double(modelFolders); %windows sorts 1 then 10 then 100 then 2 then 20, etc.

%Skip 0 time and rampTime
[~,sortID] = sort(sortModels);
skipRampTime = find(sortModels(sortID) > 20,1,'first');
sortID = sortID(skipRampTime:end);

%Model timesteps
time = modelTimeStep:modelTimeStep:(modelDuration-rampTime);
dataLength = length(sortID);

%Data bins
Zbins = -2:0.1:-0.5;

%Load the data
for i = 1:dataLength
    fileName = dir([modelBasePath modelFolders{sortID(i)} '\']);
    fileName = {fileName.name};
    isFile = contains(fileName,'epsilon');
    fid = fopen([modelBasePath modelFolders{sortID(i)} '\' fileName{isFile}]);
    file = textscan(fid,'%n%n%n%n','headerlines',2);
    x = file{1};y = file{2};z = file{3};epsilon = file{4};
    %some epsilons are negative (not sure why?). Remove these.
    epsilon(epsilon<0) = 1e-10;
    data.time(:,i) = time(i);
    data.x(:,i) = x;
    data.y(:,i) = y;
    data.z(:,i) = z;
    data.epsilon(:,i) = epsilon;
    %bin the epsilon data by the z axis
    [~,~,loc] = histcounts(z,Zbins);
    lz = find(loc>0);
    data.zBin(:,i) = accumarray(loc(lz),z(lz))./accumarray(loc(lz),1);
    data.epsilonBin(:,i) = accumarray(loc(lz),epsilon(lz))./accumarray(loc(lz),1);
    clear x y z epsilon
    fclose(fid);
end

%Plot the data
ff(1) = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[800 350   850   500]);
sp(1) = subplot(121);
p(1) = imagesc(data.time,Zbins,data.epsilonBin);
cb = colorbar;
sp(2) = subplot(122);
p(2) = plot(nanmean(data.epsilonBin,2),nanmean(data.zBin,2),...
    'linewidth',1.5,'color','k');

set(sp(1),'ydir','normal','ylim',[-2 -0.5],...
    'xlim',[0 100])
set(sp(2),'ydir','normal','ylim',[-2 -0.5],...
    'position',[0.54 0.325 0.32 0.6])
set(cb,'location','southoutside','linewidth',1.5,...
    'tickdir','out')
ylabel(sp(1),'\bf\itDepth [m]')
xlabel(sp(1),'\bf\itModel Time [s]')
xlabel(cb,'\bf\it\epsilon [W/kg]')
xlabel(sp(2),'\bf\it\epsilon [W/kg]')
suptitle(['\bf\it' modelFigName])
prettyfigures('text',11,'labels',13)
if saveFigs == 1
    export_fig(ff(1),[figPath modelName{:} '_rawTurbData_V1'],'-png','-r600','-nocrop')
end
if saveData == 1
    dataPath = 'c:\Users\user\Documents\Models\DataAnalysis\Mesh_scale_vs_turbulence\';
    save([dataPath modelName{:} '_TurbData_V1.mat'],'data','-mat','-v7.3')
end

