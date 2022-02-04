%Extract modeled turbulence values from the 'Turbulence Refinement Models'.
%This script iteratively loads the data from each models (user specified)
%into a data structure and optionally plots the results as a time series
%(image). Results can also be optionally saved to disk.
%
% Updates:
% 03/31/21: added capability to plot multiple models on the same figure
%
% This is version 2 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
modelNames = {'MKK_CombinedHR_level1';'MKK_CombinedHR_level2';'MKK_CombinedHR_level3';'MKK_CombinedHR_level4';'MKK_CombinedHR_level5'};
modelLegend = {'Level1 [0.3 - 0.15 m]';'Level2 [0.3 - 0.075 m]';'Level3 [0.3 - 0.0375 m]';...
    'Level4 [0.3 - 0.0187 m]';'Level5 [0.3 - 0.01 m]'};
%'Level7 refinement model [0.3 - 0.002 m]';
RBRbursts = [151 158]; %RBR bursts to analyze (range: 141 - 162)
rampTime = 20; %timeseries will be adjusted by this number of seconds to skip data during 'rampTime'
modelTimeStep = 0.125; %model timestep in seconds
modelDuration = 120; %length of timeseries in seconds
saveFigs = 1; %save figures to disk; 0 is off, 1 is on
saveData = 0; %save out profile data to disk; 0 is off, 1 is on

%% Plot the model data
ff(1) = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[800 350   550   500]);
hold on;
cc = brewermap(length(modelNames),'Blues');
p = zeros(length(modelNames),1);
for j = 1:length(modelNames)
    if length(RBRbursts) == 1
        burstTxt = sprintf('%d',RBRbursts(1));
    else
        burstTxt = sprintf('%d-%d',RBRbursts(1),RBRbursts(2));
    end
    modelBasePath = ['c:\Users\user\Documents\Models\MKK_Combined\2D\Burst' burstTxt '\TurbulenceRefinementModels\' modelNames{j} '\MKKmodel\postProcessing\bathySample\surface\'];
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
    Zbins = -2:0.025:-0.5;
    
    %Load the data
    for i = 1:dataLength
        fileName = dir([modelBasePath modelFolders{sortID(i)} '\']);
        fileName = {fileName.name};
        isFile = contains(fileName,'U_sampled');
        fid = fopen([modelBasePath modelFolders{sortID(i)} '\' fileName{isFile}]);
        file = textscan(fid,'%n%n%n%n%n%n','headerlines',2);
        x = file{1};y = file{2};z = file{3};Ux = file{4};Uy = file{4};Uz = file{4};
        %TEST: refine sampling area in x-direction
        cropX = find(x>=67.5 & x<=74.5);
        x = x(cropX);y = y(cropX);z = z(cropX);Ux = Ux(cropX);Uy = Uy(cropX);Uz = Uz(cropX);
        data.time(:,i) = time(i);
        data.x(:,i) = x;
        data.y(:,i) = y;
        data.z(:,i) = z; 
        data.Ux(:,i) = Ux;
        data.Uy(:,i) = Uy;
        data.Uz(:,i) = Uz;
        Umag = sqrt(Ux.^2);
        %bin the epsilon data by the z axis
        [~,~,loc] = histcounts(z,Zbins);
        lz = find(loc>0);
        data.zBin(:,i) = accumarray(loc(lz),z(lz))./accumarray(loc(lz),1);
        data.UmagBin(:,i) = accumarray(loc(lz),Umag(lz))./accumarray(loc(lz),1);
        clear x y z Ux Uy Uz
        fclose(fid);
    end
    p(j) = plot(smooth(nanmean(data.UmagBin,2),8),nanmean(data.zBin,2).*-1,...
        'linewidth',1.5,'color',cc(j,:));
    
    %Optionally save out data to disk
    if saveData == 1
        dataPath = 'c:\Users\user\Documents\Models\DataAnalysis\Mesh_scale_vs_turbulence\';
        save([dataPath modelNames{j} '_VelocityData_V2.mat'],'data','-mat','-v7.3')
    end
    clear data
end
set(gca,'ydir','reverse','ylim',[0.8 2])
leg = legend(p,modelLegend,'location','northwest');
ylabel('\bf\itDepth [m]')
xlabel('\bf\it$|\overline{U}|$ [m/s]','interpreter','latex')
prettyfigures('text',11,'labels',13)
if saveFigs == 1
    export_fig(ff(1),[figPath 'ModelScale_Umag_Profiles_V2'],'-png','-r600','-nocrop')
end


