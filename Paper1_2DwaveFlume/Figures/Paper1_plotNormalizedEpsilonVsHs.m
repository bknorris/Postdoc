% Create a plots of turbulence against normalized significant wave height 
% for HR vs LR models to directly compare the two.
%
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Load the data
%Define data paths
workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\Figures\';
HRmodelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\HR_Domain\';
LRmodelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\LR_Domain\';
% HRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\Models\';
% LRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Models\';

modelName = 'Models_allHR_LR.csv';

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelName]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);
xProfiles = 37:0.5:40;
cc = flipud(hex2rgb({'#032760';'#0057d2';'#57b0ff';'#73dcff'})); %A dark -> light blue colormap (based on Bath 112)
% cc = [32 80 255;134 217 255;255 196 0;213 0 0]./255; %A red orange light blue blue colormap (based on GMT Panopoly)
markers = {'^';'s';'o';'d'};
lineStyle = {'-';'--';'-.';':'};

%% Load processed model data -- HR
HRmodelFree = dir([HRmodelDir '*freeSurf_V2*']);
HRmodelFree = HRmodelFree(~ismember({HRmodelFree.name},{'.','..'}));
HRmodelFree = {HRmodelFree.name};
HRmodelEps = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V5.mat']);

%Loop through model scenarios and load the data into a structure
HR.WGs = [31 65 66 67 68 69 70 71 72 73 74 75 76];
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    scenarioNumber = regexp(HRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    
    %Load the free surface data
    FS = load([HRmodelDir HRmodelFree{modelID}]);
    fn = fieldnames(FS);
    Tp = modelInfo.Tp(i);
    h = modelInfo.h(i);
    for j = 1:length(HR.WGs)
        Gauges = find(FS.waveGauges == HR.WGs(j));
        [WG2,f] = pwelch(FS.(fn{Gauges}).eta,[],[],[],8);
        HR.Hs(i,j) = 4*sqrt(sum(WG2.*mean(diff(f))));
    end
    HR.Tp(i) = Tp;
    HR.h(i) = h;
    
    %Now do the epsilon data
    patchLo = 35;
    patchHi = 41;
    zPatch = [-0.9 -1.9 -2.9 -3.9];
    modelID = find(strcmp(whichModel,string(HRmodelEps.scenarioNumber)));
    eps1 = squeeze(HRmodelEps.eps(:,:,modelID));
    Uw1 = squeeze(HRmodelEps.Uw(:,:,modelID));
    xBins = squeeze(HRmodelEps.xBins(:,:,modelID));
    zBins = squeeze(HRmodelEps.zBins(:,:,modelID));
    [~,xLo] = min(abs(xBins-patchLo));
    [~,xHi] = min(abs(xBins-patchHi));
    if HRmodelEps.h(modelID) == 1
        [~,zHi] = min(abs(zBins-zPatch(1)));
    elseif HRmodelEps.h(modelID) == 2
        [~,zHi] = min(abs(zBins-zPatch(2)));
    elseif HRmodelEps.h(modelID) == 3
        [~,zHi] = min(abs(zBins-zPatch(3)));
    elseif HRmodelEps.h(modelID) == 4
        [~,zHi] = min(abs(zBins-zPatch(4)));
    end
    HR.eps(i) = nanmean(nanmean(eps1(xLo:xHi,:)));
    HR.Uw(i) = nanmean(nanmean(Uw1(xLo:xHi,:).^3));
end
sp(1) = subplot(121);
hModel = unique(HR.h);
TpModel = unique(HR.Tp);
xs = repmat(linspace(-0.6,0.6,4)',1,5);
for i = 1:length(hModel)
    hold on
    %Organize plots by water depth!
    modelID = find(HR.h == hModel(i));
    HsNorm = mean(HR.Hs(modelID,2:end)./HR.Hs(modelID,1),2);
    y = HR.eps(modelID)';
    x = HsNorm;
    
    plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);

%     edges = linspace(0,24,6);
%     bins = discretize(x,edges);
%     [meanVal,maxVal,stDev] = grpstats(y,bins,{@mean, @max, @std});
% 
%     errorbar(TpModel+xs(i,:),meanVal,stDev,...
%         sprintf('%s%s',markers{i},lineStyle{i}),'color',cc(i,:),...
%         'markerfacecolor',cc(i,:),'markersize',6,'linewidth',1.2,'capsize',0);
end

