% Create a plot of total mean TKE dissipation divided by X (patch width)
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

modelName = 'Models_allHR_LR.csv';

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelName]);
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
% %
% %Comment out when not using!
% for i = 3%1:62
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
cc = [0 57 255;0 219 218;243 245 4;219 167 0;255 3 3]./255;
markers = {'o';'s';'d';'^';'p'};
lineStyle = {'-';'--';'-.';':';'-'};

disp('Plotting HR models') %Loop through model scenarios and plot the data
hModel = unique(HR.h);
patchLo = 35;
patchHi = 46;
zPatch = [-0.9 -1.9 -2.9 -3.9];
sp1 = zeros(length(hModel),1);
spIDs = [1 3 5 7];
for i = 1:length(hModel)
    sp1(i) = subplot(4,2,spIDs(i));
    %Plot elevation of coral canopy as dashed black line
    plot(linspace(10E-6,10E-1,10),ones(1,10)*zPatch(i),'--k','linewidth',1.5)
    hold on
    
    %Organize plots by water depth!
    modelID = find(HR.h == hModel(i) & HR.Hs == 0.4);
    
    %Calculate depth profiles
    for j = 1:length(modelID)
        eps = squeeze(HR.eps(:,:,modelID(j)));
        Uw = squeeze(HR.Uw(:,:,modelID(j)));
        xBins = squeeze(HR.xBins(:,:,modelID(j)));
        zBins = squeeze(HR.zBins(:,:,modelID(j)));
        [~,xLo] = min(abs(xBins-patchLo));
        [~,xHi] = min(abs(xBins-patchHi));
        Profile = nanmean(eps(xLo:xHi,:))./nanmean(Uw(xLo:xHi,:).^3);
        z = zBins;
        x = Profile;
        plot(x,z,lineStyle{j},...
            'color',cc(j,:),'linewidth',1);
        plot(x(1:7:end),z(1:7:end),markers{j},...
            'color',cc(j,:),'markerfacecolor','w',...
            'markersize',3,'linewidth',1)
    end
    text(0.1,0.1,sprintf('h = %0.0f m',hModel(i)),'units','normalized');
end

disp('Plotting LR models') %Loop through model scenarios and plot the data
hModel = unique(LR.h);
patchLo = 35;
patchHi = 44;
zPatch = [-0.92 -1.92 -2.92 -3.92];
sp2 = zeros(length(hModel),1);
spIDs = [2 4 6 8];
for i = 1:length(hModel)
    sp2(i) = subplot(length(hModel),2,spIDs(i));
    %Plot elevation of coral canopy as dashed black line
    plot(linspace(10E-6,10E-1,10),ones(1,10)*zPatch(i),'--k','linewidth',1.5)
    hold on
    
    %Organize plots by water depth!
    modelID = find(LR.h == hModel(i) & LR.Hs == 0.4);
    
    %Calculate depth profiles
    for j = 1:length(modelID)
        eps = squeeze(LR.eps(:,:,modelID(j)));
        Uw = squeeze(LR.Uw(:,:,modelID(j)));
        xBins = squeeze(LR.xBins(:,:,modelID(j)));
        zBins = squeeze(LR.zBins(:,:,modelID(j)));
        [~,xLo] = min(abs(xBins-patchLo));
        [~,xHi] = min(abs(xBins-patchHi));
        Profile = nanmean(eps(xLo:xHi,:))./nanmean(Uw(xLo:xHi,:).^3);
        z = zBins;
        x = Profile;
        plot(x,z,lineStyle{j},...
            'color',cc(j,:),'linewidth',1);hold on;
        plot(x(1:7:end),z(1:7:end),markers{j},...
            'color',cc(j,:),'markerfacecolor','w',...
            'markersize',3,'linewidth',1)
    end
    text(0.1,0.1,sprintf('h = %0.0f m',hModel(i)),'units','normalized');
end
% 
wavePeriods = unique(LR.Tp);
%Make some dummy data for the legends
xs = repmat(linspace(1E5,1E6,2),length(wavePeriods),1);
ys = repmat(linspace(min(x),max(x),2),length(wavePeriods),1);
pp3 = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    pp3(i) = plot(xs(i,:),ys(i,:),sprintf('%s%s',lineStyle{i},markers{i}),...
        'markerfacecolor','w','color',cc(i,:),'linewidth',1);hold on;
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end
leg = legend(pp3,legText,'location','southeast','box','off','interpreter','latex');
 
%Global adjustments
set(sp1,'xscale','log','xlim',[10E-6 5E-1])
set(sp1(1),'ylim',[-1.5 0],'ytick',-1.5:0.5:0,'xticklabel',[])
set(sp1(2),'ylim',[-2.5 -1],'ytick',-2.5:0.5:-1,'xticklabel',[])
set(sp1(3),'ylim',[-3.5 -2],'ytick',-3.5:0.5:-2,'xticklabel',[])
set(sp1(4),'ylim',[-4.5 -3],'ytick',-4.5:0.5:-3)
 
set(sp2,'xscale','log','xlim',[10E-6 5E-1])
set(sp2(1),'ylim',[-1.5 0],'ytick',-1.5:0.5:0,'yticklabel',[],'xticklabel',[])
set(sp2(2),'ylim',[-2.5 -1],'ytick',-2.5:0.5:-1,'yticklabel',[],'xticklabel',[])
set(sp2(3),'ylim',[-3.5 -2],'ytick',-3.5:0.5:-2,'yticklabel',[],'xticklabel',[])
set(sp2(4),'ylim',[-4.5 -3],'ytick',-4.5:0.5:-3,'yticklabel',[])
 
%Positoning
set(sp1(1),'position',[0.1 0.76 0.34 0.2])
set(sp2(1),'position',[0.48 0.76 0.34 0.2])
set(sp1(2),'position',[0.1 0.54 0.34 0.2])
set(sp2(2),'position',[0.48 0.54 0.34 0.2])
set(sp1(3),'position',[0.1 0.32 0.34 0.2])
set(sp2(3),'position',[0.48 0.32 0.34 0.2])
set(sp1(4),'position',[0.1 0.1 0.34 0.2])
set(sp2(4),'position',[0.48 0.1 0.34 0.2])
set(leg,'position',[0.82 0.81 0.15 0.1])

%Labeling
xlabel(sp1(4),'$\langle \overline{\varepsilon}\rangle/\langle{U_w}^3 \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
xlabel(sp2(4),'$\langle \overline{\varepsilon}\rangle/\langle{U_w}^3 \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
ylabel(sp1(1),'$z$ (m)','interpreter','latex')
ylabel(sp1(2),'$z$ (m)','interpreter','latex')
ylabel(sp1(3),'$z$ (m)','interpreter','latex')
ylabel(sp1(4),'$z$ (m)','interpreter','latex')
title(sp1(1),'HR Models')
title(sp2(1),'LR Models')

prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.008 0.008])
 
set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% print(ff,[workingDir modelName(1:end-4) '_NormepsProfiles'],'-dpdf','-r0')
