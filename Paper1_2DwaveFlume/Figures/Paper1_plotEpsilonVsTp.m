% Create a plots of turbulence against wave period
% for all models
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

%Load processed model data
HR = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V5.mat']);
LR = load([LRmodelDir 'LR_postProcess_epsilon_allHR_LR_V5.mat']);

ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);
xProfiles = 37:0.5:40;
cc = flipud(hex2rgb({'#032760';'#0057d2';'#57b0ff';'#73dcff'})); %A dark -> light blue colormap (based on Bath 112)
% cc = [32 80 255;134 217 255;255 196 0;213 0 0]./255; %A red orange light blue blue colormap (based on GMT Panopoly)
markers = {'^';'s';'o';'d'};
lineStyle = {'-';'--';'-.';':'};

%% HR
sp(1) = subplot(121);
hModel = unique(HR.h);
TpModel = unique(HR.Tp);
xs = repmat(linspace(-0.6,0.6,4)',1,5);
for i = 1:length(hModel)
    hold on
    patchLo = 35;
    patchHi = 41;
    zPatch = [-0.9 -1.9 -2.9 -3.9];
    
    %Organize plots by water depth!
    modelID = find(HR.h == hModel(i));
    epsNorm = zeros(length(modelID),1);
    for j = 1:length(modelID)
        eps = squeeze(HR.eps(:,:,modelID(j)));
        Uw = squeeze(HR.Umag(:,:,modelID(j)));
        xBins = squeeze(HR.xBins(:,:,modelID(j)));
        zBins = squeeze(HR.zBins(:,:,modelID(j)));
        [~,xLo] = min(abs(xBins-patchLo));
        [~,xHi] = min(abs(xBins-patchHi));
        if HR.h(modelID(j)) == 1
            [~,zHi] = min(abs(zBins-zPatch(1)));
        elseif HR.h(modelID(j)) == 2
            [~,zHi] = min(abs(zBins-zPatch(2)));
        elseif HR.h(modelID(j)) == 3
            [~,zHi] = min(abs(zBins-zPatch(3)));
        elseif HR.h(modelID(j)) == 4
            [~,zHi] = min(abs(zBins-zPatch(4)));
        end
        epsNorm(j) = nanmean(nanmean(eps(xLo:xHi,:)./(Uw(xLo:xHi,:).^3)));
    end
    x = HR.Tp(modelID);
    y = epsNorm;

    edges = linspace(0,24,6);
    bins = discretize(x,edges);
    [meanVal,maxVal,stDev] = grpstats(y,bins,{@mean, @max, @std});

    errorbar(TpModel+xs(i,:),meanVal,stDev,...
        sprintf('%s%s',markers{i},lineStyle{i}),'color',cc(i,:),...
        'markerfacecolor',cc(i,:),'markersize',6,'linewidth',1.2,'capsize',0);
end

%% LR
sp(2) = subplot(122);
hModel = unique(LR.h);
TpModel = unique(LR.Tp);
xs = repmat(linspace(-0.6,0.6,4)',1,5);
for i = 1:length(hModel)
    hold on
    patchLo = 40.5;
    patchHi = 44;
    zPatch = [-0.92 -1.92 -2.92 -3.92];
    
    %Organize plots by water depth!
    modelID = find(LR.h == hModel(i));
    epsNorm = zeros(length(modelID),1);
    for j = 1:length(modelID)
        eps = squeeze(LR.eps(:,:,modelID(j)));
        Uw = squeeze(LR.Umag(:,:,modelID(j)));
        xBins = squeeze(LR.xBins(:,:,modelID(j)));
        zBins = squeeze(LR.zBins(:,:,modelID(j)));
        [~,xLo] = min(abs(xBins-patchLo));
        [~,xHi] = min(abs(xBins-patchHi));
        if LR.h(modelID(j)) == 1
            [~,zHi] = min(abs(zBins-zPatch(1)));
        elseif LR.h(modelID(j)) == 2
            [~,zHi] = min(abs(zBins-zPatch(2)));
        elseif LR.h(modelID(j)) == 3
            [~,zHi] = min(abs(zBins-zPatch(3)));
        elseif LR.h(modelID(j)) == 4
            [~,zHi] = min(abs(zBins-zPatch(4)));
        end
        epsNorm(j) = nanmean(nanmean(eps(xLo:xHi,:)./(Uw(xLo:xHi,:).^3)));
    end
    x = LR.Tp(modelID);
    y = epsNorm;

    edges = linspace(0,24,6);
    bins = discretize(x,edges);
    [meanVal,maxVal,stDev] = grpstats(y,bins,{@mean, @max, @std});

    errorbar(TpModel+xs(i,:),meanVal,stDev,...
        sprintf('%s%s',markers{i},lineStyle{i}),'color',cc(i,:),...
        'markerfacecolor',cc(i,:),'markersize',6,'linewidth',1.2,'capsize',0);
end

% %Make some dummy data for the legends
% xs = repmat(linspace(1E5,1E6,2),length(hModel),1);
% ys = repmat(linspace(min(x),max(x),2),length(hModel),1);
% pp = zeros(length(hModel),1);
% legText = cell(length(hModel),1);
% for i = 1:length(hModel)
%     pp(i) = plot(xs(i,:),ys(i,:),sprintf('%s%s',lineStyle{i},markers{i}),...
%         'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',6,'linewidth',1.2);hold on;
%     legText{i} = sprintf('$h = %0.0f$ m',hModel(i));
% end
% leg = legend(pp,legText,'location','southeast','box','off','interpreter','latex');
% 
%Global adjustments
set(sp(1),'xtick',4:4:20,'xlim',[2 22],...
    'ylim',[0 1.5])
set(sp(2),'xtick',4:4:20,'xlim',[2 22],...
    'ylim',[0 1.5],...
    'yticklabel',[])
% 
%Positoning
set(sp(1),'position',[0.1 0.12 0.38 0.8])
set(sp(2),'position',[0.52 0.12 0.38 0.8])
% set(leg,'position',[0.78 0.18 0.1 0.1])
% 
% %Labeling
% ylabel(sp(1),'$\langle H_s\rangle / {H_s}_0 \   \mathrm{(-)}$','interpreter','latex')
% xlabel(sp(1),'$T_p \ \mathrm{(s)}$','interpreter','latex')
% xlabel(sp(2),'$T_p \ \mathrm{(s)}$','interpreter','latex')
% title(sp(1),'HR Models')
% title(sp(2),'LR Models')
%  
% prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.008 0.008])
%  
% set(ff,'units','inches','renderer','painters');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% % print(ff,[workingDir modelName(1:end-4) '_HsVsTp_V2alt2'],'-dpdf','-r0')

