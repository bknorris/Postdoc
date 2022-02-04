% Create a color plot of time-averaged turbulence for the five depth models
% (0.5, 1, 2, 3, 4 m). Script currently plots both HR and LR models in two
% columns. 
%
% Updates:
% 06/03/21: Load the free surface time series and normalize turbulence by
% incoming wave energy dissipation (WG1 - WG2). 
%
% This is version 2 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%Define data paths
workingDir = 'c:\Users\user\Documents\Models\Figures\Paper1_2DwaveFlume\';
HRmodelDir = 'e:\BKN-FIELD\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\DataAnalysis\';
LRmodelDir = 'e:\BKN-FIELD\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\DataAnalysis\';
HRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\Models\';
LRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Models\';

modelInfo = 'Models_Hs0_4m_Tp_12_ramptime8.csv';
modelDepthScenarios = {'0.5 m';'1 m';'2 m';'3 m';'4 m'};

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Get modelNames from the model directory
HRmodelRaw = dir([HRmodelDir '*rawData*']);
HRmodelRaw = HRmodelRaw(~ismember({HRmodelRaw.name},{'.','..'}));
HRmodelRaw = {HRmodelRaw.name};

HRmodelFree = dir([HRmodelDir '*freeSurf*']);
HRmodelFree = HRmodelFree(~ismember({HRmodelFree.name},{'.','..'}));
HRmodelFree = {HRmodelFree.name};

HRmodelBathy = dir([HRbathyDir '*_bathy_0.stl']);
HRmodelBathy = HRmodelBathy(~ismember({HRmodelBathy.name},{'.','..'}));
HRmodelBathy = {HRmodelBathy.name};

LRmodelRaw = dir([LRmodelDir '*rawData_V1*']);
LRmodelRaw = LRmodelRaw(~ismember({LRmodelRaw.name},{'.','..'}));
LRmodelRaw = {LRmodelRaw.name};

LRmodelFree = dir([LRmodelDir '*freeSurf_V1*']);
LRmodelFree = LRmodelFree(~ismember({LRmodelFree.name},{'.','..'}));
LRmodelFree = {LRmodelFree.name};

LRmodelBathy = dir([LRbathyDir '*_bathy_0.stl']);
LRmodelBathy = LRmodelBathy(~ismember({LRmodelBathy.name},{'.','..'}));
LRmodelBathy = {LRmodelBathy.name};

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   850   700]);
sp1 = zeros(1,5);
HRsubplots = [1 3 5 7 9];
disp('Loading HR models')
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(HRmodelRaw','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),HRmodelRaw{modelID})
    load([HRmodelDir HRmodelRaw{modelID}]);
    HR = data;clear data
    load([HRmodelDir HRmodelFree{modelID}]);
    FS = data;clear data
    
    %Calculate wave energy dissipation @ WG1 and WG2 
    [WG1,~] = pwelch(FS.WG1.eta,[],[],[],8);
    [WG2,F] = pwelch(FS.WG2.eta,[],[],[],8);
    
    k=wavek(F,modelInfo.h(i));
    kh = k*modelInfo.h(i);
    c = sqrt(9.81.*tanh(kh))./sqrt(k);
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c; 
    F1 = trapz(1025*9.81*WG1(2:end).*Cg(2:end));
    F2 = trapz(1025*9.81*WG2(2:end).*Cg(2:end));
    D = (F1-F2)/(63-31); %wave energy dissipation between WG1 and WG2 (Huang et al. 2012)
    
    %Average turbulence in time
    epsAvg = mean(HR.epsilon.epsilon,2);
    x = HR.epsilon.x(:,1)-30;
    z = HR.epsilon.z(:,1);
    
    %Normalize turbulence by incoming wave energy dissipation
    epsNorm = epsAvg./D;

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
    eBins = linspace(min(epsNorm),max(epsNorm),nBins);
    
    ix = discretize(x, xBins);
    iz = discretize(z, zBins);
    idx = sub2ind([nBins nBins], ix, iz);
    epsg = accumarray(idx, epsNorm, [nBins*nBins 1], @(x) mean(x), NaN);
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
    plot(bathyX-30,bathyY,'-k','linewidth',2)

    set(gca,'ydir','normal')
    cc = brewermap(nBins,'*RdYlBu');
    colormap(cc)
    clear HR
    caxis([1e-9 1e-5])
end

sp2 = zeros(1,5);
LRsubplots = [2 4 6 8 10];
disp('Loading LR models')
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(LRmodelRaw','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),LRmodelRaw{modelID})
    load([LRmodelDir LRmodelRaw{modelID}]);
    LR = data;clear data
    load([LRmodelDir LRmodelFree{modelID}]);
    FS = data;clear data
    
    %Calculate wave energy dissipation @ WG1 and WG2 
    [WG1,~] = pwelch(FS.WG1.eta,[],[],[],8);
    [WG2,F] = pwelch(FS.WG2.eta,[],[],[],8);
    
    k=wavek(F,modelInfo.h(i));
    kh = k*modelInfo.h(i);
    c = sqrt(9.81.*tanh(kh))./sqrt(k);
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = trapz(n(2:end))*trapz(c(2:end)); 
    F1 = trapz(WG1(2:end))*Cg*(9.81)*1025;
    F2 = trapz(WG2(2:end))*Cg*(9.81)*1025;
    D = (F1-F2)/(63-31); %wave energy dissipation between WG1 and WG2 (Huang et al. 2012)
    
    %Average turbulence in time
    epsAvg = mean(LR.epsilon.epsilon,2);
    x = LR.epsilon.x(:,1)-30;
    z = LR.epsilon.z(:,1);
    
    %Normalize turbulence by incoming wave energy dissipation
    epsNorm = epsAvg./D;

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
    eBins = linspace(min(epsNorm),max(epsNorm),nBins);
    
    ix = discretize(x, xBins);
    iz = discretize(z, zBins);
    idx = sub2ind([nBins nBins], ix, iz);
    epsg = accumarray(idx, epsNorm, [nBins*nBins 1], @(x) mean(x), NaN);
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
    caxis([1e-9 1e-5])
end

%Create legend
cb1 = colorbar(sp1(5),'location','southoutside');
set(cb1,'position',[0.1 0.08 0.38 0.02],...
    'tickdir','out',...
    'ticklength',0.01,...
    'linewidth',1.5);
cb2 = colorbar(sp2(5),'location','southoutside');
set(cb2,'position',[0.55 0.08 0.38 0.02],...
    'tickdir','out',...
    'ticklength',0.01,...
    'linewidth',1.5);

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
xlabel(cb1,'$\overline{\epsilon}/\overline{\epsilon_0} \quad (-)$','interpreter','latex')
xlabel(cb2,'$\overline{\epsilon}/\overline{\epsilon_0} \quad (-)$','interpreter','latex')

%
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.01 0.01])

%Save figure
% export_fig(ff,[workingDir  'HR_LR_Tp_12s_TimeAvgd_Turbulence_V3'],'-png','-r600','-nocrop')