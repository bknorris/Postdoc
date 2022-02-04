% Create a color plot of time-averaged turbulence for the five depth models
% (0.5, 1, 2, 3, 4 m). Script currently plots both HR and LR models in two
% columns. 
%
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%Define data paths
workingDir = 'c:\Users\user\Documents\Models\Figures\Paper1_2DwaveFlume\';
HRmodelDir = 'f:\USGS\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\DataAnalysis\';
LRmodelDir = 'f:\USGS\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\DataAnalysis\';
HRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\Models\';
LRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Models\';

modelInfo = 'Models_Hs_0_4m_Tp_4-12.csv';
modelDepthScenarios = {'0.5 m';'1 m';'2 m';'3 m';'4 m'};

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Get modelNames from the model directory
HRmodelScenarios = dir([HRmodelDir '*rawData*']);
HRmodelScenarios = HRmodelScenarios(~ismember({HRmodelScenarios.name},{'.','..'}));
HRmodelScenarios = {HRmodelScenarios.name};

HRmodelBathy = dir([HRbathyDir '*_bathy_0.stl']);
HRmodelBathy = HRmodelBathy(~ismember({HRmodelBathy.name},{'.','..'}));
HRmodelBathy = {HRmodelBathy.name};


LRmodelScenarios = dir([LRmodelDir '*rawData*']);
LRmodelScenarios = LRmodelScenarios(~ismember({LRmodelScenarios.name},{'.','..'}));
LRmodelScenarios = {LRmodelScenarios.name};

LRmodelBathy = dir([LRbathyDir '*_bathy_0.stl']);
LRmodelBathy = LRmodelBathy(~ismember({LRmodelBathy.name},{'.','..'}));
LRmodelBathy = {LRmodelBathy.name};

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   850   700]);
sp1 = zeros(1,10);
HRsubplots = [1 3 5 7 9];
disp('Loading HR models')
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(HRmodelScenarios','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),HRmodelScenarios{modelID})
    load([HRmodelDir HRmodelScenarios{modelID}]);
    HR = data;clear data
    
    %Average turbulence in time
    epsAvg = mean(HR.epsilon.epsilon,2);
    x = HR.epsilon.x(:,1);
    z = HR.epsilon.z(:,1);

%     bins = linspace(min(epsAvg),max(epsAvg),1000);
%     [N,edges,histIDs] = histcounts(epsAvg,bins);
%     cc = brewermap(length(N),'*RdYlBu');
%     for j = 1:length(N) 
%         eachBin = find(histIDs == j);
%         plot3(x(eachBin),z(eachBin),epsAvg(eachBin),'.','color',cc(j,:))
%         hold on
%     end
%     view(0,90)
    
    nBins = 100;
    xBins = linspace(min(x),max(x),nBins);
    zBins = linspace(min(z),max(z),nBins);
    eBins = linspace(min(epsAvg),max(epsAvg),nBins);
    
    ix = discretize(x, xBins);
    iz = discretize(z, zBins);
    idx = sub2ind([nBins nBins], ix, iz);
    epsg = accumarray(idx, epsAvg, [nBins*nBins 1], @(x) mean(x), NaN);
    epsg = reshape(epsg, nBins, nBins);
    
    sp1(i) = subplot(5,2,HRsubplots(i));
    %First plot the data!
    p = pcolor(xBins,zBins,epsg');
    shading interp
    hold on
    %Now plot the bathy as a black line
    data = stlread([HRbathyDir HRmodelBathy{i}]);
    [bathyX,sortID] = sort(data.Points(:,1));
    bathyY = data.Points(sortID,3);
    plot(bathyX,bathyY,'-k','linewidth',2)

    set(gca,'ydir','normal')
    cc = brewermap(nBins,'*RdYlBu');
    colormap(cc)
    clear HR
    caxis([1e-6 1e-4])
end

sp2 = zeros(1,10);
LRsubplots = [2 4 6 8 10];
disp('Loading LR models')
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(LRmodelScenarios','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),LRmodelScenarios{modelID})
    load([LRmodelDir LRmodelScenarios{modelID}]);
    LR = data;clear data
    
    %Average turbulence in time
    epsAvg = mean(LR.epsilon.epsilon,2);
    x = LR.epsilon.x(:,1)-30;
    z = LR.epsilon.z(:,1);

%     bins = linspace(min(epsAvg),max(epsAvg),1000);
%     [N,edges,histIDs] = histcounts(epsAvg,bins);
%     cc = brewermap(length(N),'*RdYlBu');
%     for j = 1:length(N) 
%         eachBin = find(histIDs == j);
%         plot3(x(eachBin),z(eachBin),epsAvg(eachBin),'.','color',cc(j,:))
%         hold on
%     end
%     view(0,90)
    
    nBins = 100;
    xBins = linspace(min(x),max(x),nBins);
    zBins = linspace(min(z),max(z),nBins);
    eBins = linspace(min(epsAvg),max(epsAvg),nBins);
    
    ix = discretize(x, xBins);
    iz = discretize(z, zBins);
    idx = sub2ind([nBins nBins], ix, iz);
    epsg = accumarray(idx, epsAvg, [nBins*nBins 1], @(x) mean(x), NaN);
    epsg = reshape(epsg, nBins, nBins);
    
    sp2(i) = subplot(5,2,LRsubplots(i));
    %First plot the data!
    p = pcolor(xBins,zBins,epsg');
    shading interp
    hold on
    %Now plot the bathy as a black line
    data = stlread([LRbathyDir LRmodelBathy{i}]);
    [bathyX,sortID] = sort(data.Points(:,1));
    bathyY = data.Points(sortID,3);
    plot(bathyX-30,bathyY,'-k','linewidth',2)

    set(gca,'ydir','normal')
    cc = brewermap(nBins,'*RdYlBu');
    colormap(cc)
    clear LR
    caxis([1e-6 1e-4])
end

%Create legend
cb = colorbar('location','southoutside');
set(cb,'position',[0.16 0.08 0.7 0.02],...
    'tickdir','out',...
    'ticklength',0.01,...
    'linewidth',1.5)

% %Global Adjustments
set([sp1(1) sp2(1)],'ylim',[-1.1 -0.2],'ytick',-1:0.4:-0.2)
set([sp1(2) sp2(2)],'ylim',[-1.6 -0.4])
set([sp1(3) sp2(3)],'ylim',[-2.6 -1])
set([sp1(4) sp2(4)],'ylim',[-3.6 -2])
set([sp1(5) sp2(5)],'ylim',[-4.6 -3])

set(sp1(1),'position',[0.1 0.8 0.38 0.12],'xticklabel',[])
set(sp1(2),'position',[0.1 0.65 0.38 0.12],'xticklabel',[])
set(sp1(3),'position',[0.1 0.5 0.38 0.12],'xticklabel',[])
set(sp1(4),'position',[0.1 0.35 0.38 0.12],'xticklabel',[])
set(sp1(5),'position',[0.1 0.2 0.38 0.12])

set(sp2(1),'position',[0.55 0.8 0.38 0.12],'xticklabel',[],'yticklabel',[])
set(sp2(2),'position',[0.55 0.65 0.38 0.12],'xticklabel',[],'yticklabel',[])
set(sp2(3),'position',[0.55 0.5 0.38 0.12],'xticklabel',[],'yticklabel',[])
set(sp2(4),'position',[0.55 0.35 0.38 0.12],'xticklabel',[],'yticklabel',[])
set(sp2(5),'position',[0.55 0.2 0.38 0.12],'yticklabel',[])

% %Labeling
title(sp1(1),'HR Models')
title(sp2(1),'LR Models')
ylabel(sp1(1),'z (m)')
ylabel(sp1(2),'z (m)')
ylabel(sp1(3),'z (m)')
ylabel(sp1(4),'z (m)')
ylabel(sp1(5),'z (m)')
xlabel(sp1(5),'Across-shore Distance (m)')
xlabel(sp2(5),'Across-shore Distance (m)')
xlabel(cb,'$\overline{\epsilon} \quad (m^2/s^3)$','interpreter','latex')
%
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.01 0.01])

%Save figure
export_fig(ff,[workingDir  'HR_LR_Tp_4s_TimeAvgd_Turbulence'],'-png','-r600','-nocrop')