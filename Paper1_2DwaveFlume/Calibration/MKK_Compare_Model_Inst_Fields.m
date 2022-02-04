%Compare model derived quantities (e.g., U, epsilon) with AQDP measurements.
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
modelBasePath = 'c:\Users\user\Documents\Models\MKK_Combined\2D\MKK_CombinedHR_V15\MKKmodel\postProcessing\line\';
% modelBasePath = 'c:\Users\user\Documents\Models\MKK_Combined\2D\MKK_CombinedLR_V1\MKKmodel\postProcessing\line\';
lineName = {'line2_alpha.water_epsilon_p_rgh.xy';'line2_U.xy'};
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
insts = 'MKK18HR101aqdHR_waves.mat';
% insts = 'MKK18LR101aqdHR_waves.mat';
modelName = 'MKK_CombinedHR_V15';
% modelName = 'MKK_CombinedLR_V1';

%% Extract data from the model directory
rampTime = 128; %find values greater than rampTime
modelFolders = dir(modelBasePath);
idx = cellfun(@(x)str2double(x) > rampTime,{modelFolders.name},'UniformOutput',false);
idx = cell2mat(idx);
idx2 = find(idx);
eM = zeros(100,length(modelFolders)-3); %model epsilon
uM = zeros(100,length(modelFolders)-3); %model U
vM = zeros(100,length(modelFolders)-3); %model V
wM = zeros(100,length(modelFolders)-3); %model W
for i = 4:length(modelFolders) %skip . and ..
    fileName = dir([modelBasePath modelFolders(i).name '\']);
    fileName = {fileName.name};
    %Load epsilon data
    isFile = contains(fileName,lineName{1});
    fid = fopen([modelBasePath modelFolders(i).name '\' fileName{isFile}]);
    data = textscan(fid,'%n%n%n%n');
    freeS = find(data{2} >= 0.5);
    eM(1:length(freeS),i) = data{3}(freeS);
    fclose(fid);
    clear data
    
    %Load U data
    isFile = contains(fileName,lineName{2});
    fid = fopen([modelBasePath modelFolders(i).name '\' fileName{isFile}]);
    data = textscan(fid,'%n%n%n%n');
    uM(1:length(freeS),i) = data{2}(freeS);
    vM(1:length(freeS),i) = data{3}(freeS);
    wM(1:length(freeS),i) = data{4}(freeS);
    fclose(fid);
    clear data
end

%% Now load the instrument data
load([dataPath insts]);
whichBurst = 54; %burst to perform analysis on (hint: check MKK_ModelLogs.xlsx)
burstID = find(wave.burst == whichBurst);

%Find depths consistent between the model and the instrument
depthM = linspace(0,2.35,100); %model depths same length as xy line
sameDepth = find(depthM <= max(wave.zuv));

%Plot some simple comparisons for now
ff(1) = figure(1);
    set(ff,'PaperOrientation','landscape',...
    'position',[800 350   1250   500]);
sp(1) = subplot(141);
hold on
p(1) = plot(wave.Urot(:,burstID),wave.zuv,'linewidth',1.5,'color','b');
p(2) = plot(mean(uM(sameDepth,:),2),depthM(sameDepth),'linewidth',1.5,'color','r');
hold off
sp(2) = subplot(142);
hold on
p(1) = plot(wave.Vrot(:,burstID),wave.zuv,'linewidth',1.5,'color','b');
p(2) = plot(mean(vM(sameDepth,:),2),depthM(sameDepth),'linewidth',1.5,'color','r');
hold off
sp(3) = subplot(143);
hold on
p(1) = plot(wave.W(:,burstID),wave.zuv,'linewidth',1.5,'color','b');
p(2) = plot(mean(wM(sameDepth,:),2),depthM(sameDepth),'linewidth',1.5,'color','r');
hold off
sp(4) = subplot(144);
hold on
p(1) = plot(wave.epsilon(:,burstID),wave.zuv(1:14),'linewidth',1.5,'color','b');
p(2) = plot(mean(eM(sameDepth,:),2),depthM(sameDepth),'linewidth',1.5,'color','r');
hold off
leg = legend(p,{'AQDP';'Model'});
ylabel(sp(1),'Depth [m]')
xlabel(sp(1),'U [m/s]')
xlabel(sp(2),'V [m/s]')
xlabel(sp(3),'W [m/s]')
xlabel(sp(4),'\epsilon [m^2/s^2]')
set(sp,'ytick',0:0.1:0.7)
suptitle(sprintf('Burst no. %d',whichBurst))
prettyfigures('text',11,'labels',13)
export_fig(ff(1),[figPath modelName '_AQDP_MODEL-Fields'],'-png')

