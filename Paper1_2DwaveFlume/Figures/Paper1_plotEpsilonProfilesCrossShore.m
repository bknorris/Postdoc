% Create a plot of TKE dissipation per cross-shore distance (at fixed x)
% for each wave period model.
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
HRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\Models\';
LRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Models\';

modelInfo = 'Models_allHR_LR.csv';

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);clear data

%Load processed model data
HR = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V2.mat']);
LR = load([LRmodelDir 'LR_postProcess_epsilon_allHR_LR_V2.mat']);

HRmodelBathy = dir([HRbathyDir '*_bathy_0.stl']);
HRmodelBathy = HRmodelBathy(~ismember({HRmodelBathy.name},{'.','..'}));
HRmodelBathy = {HRmodelBathy.name};

LRmodelBathy = dir([LRbathyDir '*_bathy_0.stl']);
LRmodelBathy = LRmodelBathy(~ismember({LRmodelBathy.name},{'.','..'}));
LRmodelBathy = {LRmodelBathy.name};

%% Test plots

%Comment out when not using!
% for i = 1%1:62
%     ff = figure(1);
%     set(ff,'PaperOrientation','landscape',...
%     'position',[500 80   650   500]);
%     modelID = i;
%     fprintf('Scenario number: %0.0f\n',HR.scenarioNumber(modelID))
%     
%     xBins = squeeze(HR.xBins(:,:,modelID));
%     zBins = squeeze(HR.zBins(:,:,modelID));
%     eps = squeeze(HR.eps(:,:,modelID));
%     imagesc(xBins,zBins,eps')
%     set(gca,'ydir','normal')
%     hold on
%     %Now plot the bathy as a black line
%     hModel = HR.h(modelID);
%     hBathy = regexp(HRmodelBathy,'(?<=_)(.*?)(?=m)','match');
%     hBathy = str2double([hBathy{:}]);
%     bathyID = find(hModel == hBathy);
%     data = stlread([HRbathyDir HRmodelBathy{bathyID}]);
%     [bathyX,sortID] = sort(data.Points(:,1));
%     bathyY = data.Points(sortID,3);
%     plot(bathyX-30,bathyY,'-k','linewidth',2)
%     set(gca,'ydir','normal')
%     cc = brewermap(length(zBins),'*RdYlBu');
%     colormap(cc)
%     colorbar;
%     caxis([10E-8 10E-3])
% %     pause()
% %     close(ff)
% end
% 
% % %Comment out when not using!
% for i = 1%1:62
%     ff = figure(2);
%     set(ff,'PaperOrientation','landscape',...
%     'position',[500 80   650   500]);
%     modelID = i;
%     fprintf('Scenario number: %0.0f\n',LR.scenarioNumber(modelID))
%     
%     xBins = squeeze(LR.xBins(:,:,modelID));
%     zBins = squeeze(LR.zBins(:,:,modelID));
%     eps = squeeze(LR.eps(:,:,modelID));
%     imagesc(xBins,zBins,eps')
%     set(gca,'ydir','normal')
%     hold on
%     %Now plot the bathy as a black line
%     hModel = LR.h(modelID);
%     hBathy = regexp(LRmodelBathy,'(?<=_)(.*?)(?=m)','match');
%     hBathy = str2double([hBathy{:}]);
%     bathyID = find(hModel == hBathy);
%     data = stlread([LRbathyDir LRmodelBathy{bathyID}]);
%     [bathyX,sortID] = sort(data.Points(:,1));
%     bathyY = data.Points(sortID,3);
%     plot(bathyX-30,bathyY,'-k','linewidth',2)
%     set(gca,'ydir','normal')
%     cc = brewermap(length(zBins),'*RdYlBu');
%     colormap(cc)
%     colorbar;
%     caxis([10E-8 10E-3])
% %     pause()
% %     close(ff)
% end

%% Create Figure
ff = figure(2);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   850   1000]);
disp('Plotting HR models')
hModel = unique(HR.h);
xProfiles = 37:0.5:40;
scale = 5E-4;
sp1 = zeros(length(hModel),1);
cc = [0    0.2235    1.0000
    0.3977    0.4007    0.7874
    0.7900    0.5696    0.5736
    0.8965    0.2930    0.2930
    1.0000    0.0118    0.0118];
markers = {'o';'s';'d';'^';'p'};
%Loop through model scenarios and plot the data
spIDs = [1 3 5 7];
for i = 1:length(hModel)
    sp1(i) = subplot(4,2,spIDs(i));
    
    %Organize plots by water depth!
    modelID = find(HR.h == hModel(i) & HR.Hs == 0.4);
    
    %First plot bathymetry
    hBathy = regexp(HRmodelBathy,'(?<=_)(.*?)(?=m)','match');
    hBathy = str2double([hBathy{:}]);
    bathyID = find(hModel(i) == hBathy);
    data = stlread([HRbathyDir HRmodelBathy{bathyID}]);
    [bathyX,sortID] = sort(data.Points(:,1));
    bathyY = data.Points(sortID,3);
    plot(bathyX-30,bathyY,'-k','linewidth',2)
    hold on;
    
    %Calculate depth profiles
    for j = 1:length(modelID)
        eps = squeeze(HR.eps(:,:,modelID(j)));
        D = HR.D(modelID(j));
        xBins = squeeze(HR.xBins(:,:,modelID(j)));
        zBins = squeeze(HR.zBins(:,:,modelID(j)));
        for k = 1:length(xProfiles)
            [~,xLoc] = min(abs(xProfiles(k)-xBins));
            epsProfile = eps(xLoc,:)./D;
            z = zBins;
            x = xProfiles(k)+(epsProfile./scale);
            plot(x,z,'-',...
            'color',cc(j,:),'linewidth',1);hold on;
            plot(x(1:7:end),z(1:7:end),'o',...
            'color',cc(j,:),'markerfacecolor','w',...
            'markersize',2,'linewidth',1)
        end
    end
end
set(sp1,'xlim',[36.5 40.5])
set(sp1(1),'ylim',[-1.5 -0.5])
set(sp1(2),'ylim',[-2.5 -1.5])
set(sp1(3),'ylim',[-3.5 -2.5])
set(sp1(4),'ylim',[-4.5 -3.5])
disp('Plotting LR models')
hModel = unique(LR.h);
xProfiles = 32:0.5:35;
scale = 5E-4;
sp2 = zeros(length(hModel),1);
spIDs = [2 4 6 8];
for i = 1:length(hModel)
    sp2(i) = subplot(4,2,spIDs(i));
    
    %Organize plots by water depth!
    modelID = find(LR.h == hModel(i) & LR.Hs == 0.4);
    
    %First plot bathymetry
    hBathy = regexp(LRmodelBathy,'(?<=_)(.*?)(?=m)','match');
    hBathy = str2double([hBathy{:}]);
    bathyID = find(hModel(i) == hBathy);
    data = stlread([LRbathyDir LRmodelBathy{bathyID}]);
    [bathyX,sortID] = sort(data.Points(:,1));
    bathyY = data.Points(sortID,3);
    plot(bathyX-30,bathyY,'-k','linewidth',2)
    hold on;
    
    %Calculate depth profiles
    for j = 1:length(modelID)
        eps = squeeze(LR.eps(:,:,modelID(j)));
        D = LR.D(modelID(j));
        xBins = squeeze(LR.xBins(:,:,modelID(j)));
        zBins = squeeze(LR.zBins(:,:,modelID(j)));
        for k = 1:length(xProfiles)
            [~,xLoc] = min(abs(xProfiles(k)-xBins));
            epsProfile = eps(xLoc,:)./D;
            z = zBins;
            x = xProfiles(k)+(epsProfile./scale);
            plot(x,z,'-',...
            'color',cc(j,:),'linewidth',1);hold on;
            plot(x(1:7:end),z(1:7:end),'o',...
            'color',cc(j,:),'markerfacecolor','w',...
            'markersize',2,'linewidth',1)
        end
    end
end
set(sp2,'xlim',[31.5 35.5])
set(sp2(1),'ylim',[-1.5 -0.5])
set(sp2(2),'ylim',[-2.5 -1.5])
set(sp2(3),'ylim',[-3.5 -2.5])
set(sp2(4),'ylim',[-4.5 -3.5])

