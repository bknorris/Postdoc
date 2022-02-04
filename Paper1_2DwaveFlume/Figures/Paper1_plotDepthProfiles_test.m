% Create a plot of mean U and turbulent intensity (sqrt(k)/Umag) for
% various normalized cross-shore positions within the patches
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


%% Create Figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   950   1000]);

disp('Plotting HR models') %Loop through model scenarios and plot the data
hModel = unique(HR.h);
TpModel = unique(HR.Tp);
patchLo = 35;
patchHi = 46;
zPatch = [-0.9 -1.9 -2.9 -3.9];
sp = zeros(length(hModel),3);

%First plot pcolor maps of the Tp = 12 models for each depth
spIDs = [1 4 7 10];
for i = 1:length(hModel)
    sp(i,1) = subplot(4,3,spIDs(i));
    
    %Organize plots by water depth!
    modelID = find(HR.h == hModel(i) & HR.Hs == 0.4 & HR.Tp == 12);
    
    xBins = squeeze(HR.xBins(:,:,modelID));
    zBins = squeeze(HR.zBins(:,:,modelID));
    [~,xLo] = min(abs(xBins-patchLo));
    [~,xHi] = min(abs(xBins-patchHi));
    eps = squeeze(HR.eps(:,:,modelID));
    pcolor(xBins(xLo:xHi),zBins,eps(xLo:xHi,:)')
    shading interp
    set(gca,'ydir','normal')
    hold on
    %Now plot the bathy as a black line
    hBathy = regexp(LRmodelBathy,'(?<=_)(.*?)(?=m)','match');
    hBathy = str2double([hBathy{:}]);
    bathyID = find(hModel(i) == hBathy);
    data = stlread([HRbathyDir HRmodelBathy{bathyID}]);
    [bathyX,sortID] = sort(data.Points(:,1));
    bathyY = data.Points(sortID,3);
    plot(bathyX-30,bathyY-0.01,'-k','linewidth',1)
    set(gca,'ydir','normal')
    cc1 = brewermap(length(zBins),'Reds');
    colormap(cc1)
    caxis([10E-6 10E-3])
    if i == 4
        cb = colorbar;
        set(cb,'location','southoutside')
    end
end

%Next plot depth profiles of Uw (wave orbital velocity)
cc2 = [0 57 255;0 219 218;243 245 4;219 167 0;255 3 3]./255;
markers = {'o';'s';'d';'^';'p'};
lineStyle = {'-';'--';'-.';':';'-'};
spIDs = [2 5 8 11];
for i = 1:length(hModel)
    sp(i,2) = subplot(4,3,spIDs(i));
      
    %Plot elevation of coral canopy as dashed black line
    plot(linspace(10E-6,10E1,10),ones(1,10)*zPatch(i),'--k','linewidth',1.5)
    hold on
    
    for j = 1:length(TpModel)
        %Organize plots by water depth then by wave period!
        modelID = find(HR.h == hModel(i) & HR.Tp == TpModel(j));
        xs = zeros(100,length(modelID));
        zs = zeros(100,length(modelID));
        %Calculate depth profiles
        for k = 1:length(modelID)
            Uw = squeeze(HR.Uw(:,:,modelID(k)));
            xBins = squeeze(HR.xBins(:,:,modelID(k)));
            zBins = squeeze(HR.zBins(:,:,modelID(k)));
            [~,xLo] = min(abs(xBins-patchLo));
            [~,xHi] = min(abs(xBins-patchHi));
            Profile = nanmean(Uw(xLo:xHi,:),1);
            zs(:,k) = zBins;
            xs(:,k)  = Profile;
        end
        xsMean = nanmean(xs,2);
        xsStd = std(xs,1,2);
%         errorbar(xsMean(1:7:end),zs(1:7:end,1),xsStd(1:7:end),'horizontal',...
%             sprintf('%s',markers{j}),'color',cc2(j,:),...
%             'markerfacecolor',cc2(j,:),'markersize',5,'linewidth',1);
        plot(xsMean,zs(:,1),lineStyle{j},...
            'color',cc2(j,:),'linewidth',1);
        plot(xsMean(1:7:end),zs(1:7:end,:),markers{j},...
            'color',cc2(j,:),'markerfacecolor','w',...
            'markersize',3,'linewidth',1)
    end
end

%Lastly plot depth profiles of normalized turbulence (eps/Uw^3)
spIDs = [3 6 9 12];
for i = 1:length(hModel)
    sp(i,3) = subplot(4,3,spIDs(i));
      
    %Plot elevation of coral canopy as dashed black line
    plot(linspace(10E-6,10E1,10),ones(1,10)*zPatch(i),'--k','linewidth',1.5)
    hold on
    
    %Organize plots by water depth!
    modelID = find(HR.h == hModel(i) & HR.Hs == 0.4);
    
    for j = 1:length(TpModel)
        %Organize plots by water depth then by wave period!
        modelID = find(HR.h == hModel(i) & HR.Tp == TpModel(j));
        xs = zeros(100,length(modelID));
        zs = zeros(100,length(modelID));
        %Calculate depth profiles
        for k = 1:length(modelID)
            eps = squeeze(HR.eps(:,:,modelID(k)));
            Uw = squeeze(HR.Uw(:,:,modelID(k)));
            xBins = squeeze(HR.xBins(:,:,modelID(k)));
            zBins = squeeze(HR.zBins(:,:,modelID(k)));
            [~,xLo] = min(abs(xBins-patchLo));
            [~,xHi] = min(abs(xBins-patchHi));
            Profile = nanmean(eps(xLo:xHi,:))./nanmean(Uw(xLo:xHi,:).^3);
            zs(:,k) = zBins;
            xs(:,k)  = Profile;
        end
        xsMean = nanmean(xs,2);
        xsStd = std(xs,1,2);
%         errorbar(xsMean(1:7:end),zs(1:7:end,1),xsStd(1:7:end),'horizontal',...
%             sprintf('%s',markers{j}),'color',cc2(j,:),...
%             'markerfacecolor',cc2(j,:),'markersize',5,'linewidth',1);
        plot(xsMean,zs(:,1),lineStyle{j},...
            'color',cc2(j,:),'linewidth',1);
        plot(xsMean(1:7:end),zs(1:7:end,:),markers{j},...
            'color',cc2(j,:),'markerfacecolor','w',...
            'markersize',3,'linewidth',1)
    end
end
wavePeriods = unique(HR.Tp);
%Make some dummy data for the legends
xs = repmat(linspace(1E5,1E6,2),length(wavePeriods),1);
ys = repmat(linspace(10,100,2),length(wavePeriods),1);
pp = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    pp(i) = plot(xs(i,:),ys(i,:),sprintf('%s%s',lineStyle{i},markers{i}),...
        'markerfacecolor','w','color',cc2(i,:),'linewidth',1);hold on;
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end
leg = legend(pp,legText,'location','southeast','box','on','interpreter','latex');

%Global adjustments
set(sp(:,3),'xscale','log','xlim',[10E-6 10E0],'yticklabel',[])
set(sp(:,2),'xlim',[0 1.5],'yticklabel',[])
set(sp(1,:),'ylim',[-1.75 -0.25],'ytick',-1.5:0.5:0,'xticklabel',[])
set(sp(2,:),'ylim',[-2.75 -1.25],'ytick',-2.5:0.5:-1,'xticklabel',[])
set(sp(3,:),'ylim',[-3.75 -2.25],'ytick',-3.5:0.5:-2,'xticklabel',[])
set(sp(4,:),'ylim',[-4.75 -3.25],'ytick',-4.5:0.5:-3)
 
%Positoning
set(sp(1,1),'position',[0.1 0.78 0.25 0.2])
set(sp(1,2),'position',[0.37 0.78 0.25 0.2])
set(sp(1,3),'position',[0.64 0.78 0.25 0.2])

set(sp(2,1),'position',[0.1 0.56 0.25 0.2])
set(sp(2,2),'position',[0.37 0.56 0.25 0.2])
set(sp(2,3),'position',[0.64 0.56 0.25 0.2])

set(sp(3,1),'position',[0.1 0.34 0.25 0.2])
set(sp(3,2),'position',[0.37 0.34 0.25 0.2])
set(sp(3,3),'position',[0.64 0.34 0.25 0.2])

set(sp(4,1),'position',[0.1 0.12 0.25 0.2])
set(sp(4,2),'position',[0.37 0.12 0.25 0.2])
set(sp(4,3),'position',[0.64 0.12 0.25 0.2])


set(leg,'position',[0.86 0.84 0.1 0.1])
set(cb,'position',[0.1 0.05 0.25 0.015],'linewidth',1)

%Labeling
xlabel(cb,'$\overline{\varepsilon} \  \mathrm{(m^2/s^3)}$','interpreter','latex')
xlabel(sp(4,1),'$\mathrm{Across \ Shore \ Distance \ (m)}$','interpreter','latex')
xlabel(sp(4,2),'$\langle \overline{U_w}\rangle \  \mathrm{(m/s)}$','interpreter','latex')
xlabel(sp(4,3),'$\langle \overline{\varepsilon}\rangle/\langle{U_w}^3 \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
ylabel(sp(1,1),'$z$ (m)','interpreter','latex')
ylabel(sp(2,1),'$z$ (m)','interpreter','latex')
ylabel(sp(3,1),'$z$ (m)','interpreter','latex')
ylabel(sp(4,1),'$z$ (m)','interpreter','latex')

prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.008 0.008])
 
set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% print(ff,[workingDir modelName(1:end-4) '_HRdepthProfiles'],'-dpdf','-r0')

% %% Create Figure
% ff = figure(2);
% set(ff,'PaperOrientation','landscape',...
%     'position',[100 80   950   1000]);
% 
% disp('Plotting LR models') %Loop through model scenarios and plot the data
% hModel = unique(LR.h);
% patchLo = 35;
% patchHi = 44;
% zPatch = [-0.92 -1.92 -2.92 -3.92];
% sp = zeros(length(hModel),3);
% 
% %First plot pcolor maps of the Tp = 12 models for each depth
% spIDs = [1 4 7 10];
% for i = 1:length(hModel)
%     sp(i,1) = subplot(4,3,spIDs(i));
%     
%     %Organize plots by water depth!
%     modelID = find(LR.h == hModel(i) & LR.Hs == 0.4 & LR.Tp == 12);
%     
%     xBins = squeeze(LR.xBins(:,:,modelID));
%     zBins = squeeze(LR.zBins(:,:,modelID));
%     [~,xLo] = min(abs(xBins-patchLo));
%     [~,xHi] = min(abs(xBins-patchHi));
%     eps = squeeze(LR.eps(:,:,modelID));
%     pcolor(xBins(xLo:xHi),zBins,eps(xLo:xHi,:)')
%     shading interp
%     set(gca,'ydir','normal')
%     hold on
%     %Now plot the bathy as a black line
%     hBathy = regexp(LRmodelBathy,'(?<=_)(.*?)(?=m)','match');
%     hBathy = str2double([hBathy{:}]);
%     bathyID = find(hModel(i) == hBathy);
%     data = stlread([LRbathyDir LRmodelBathy{bathyID}]);
%     [bathyX,sortID] = sort(data.Points(:,1));
%     bathyY = data.Points(sortID,3);
%     plot(bathyX-30,bathyY-0.01,'-k','linewidth',1)
%     set(gca,'ydir','normal')
%     cc1 = brewermap(length(zBins),'Reds');
%     colormap(cc1)
%     caxis([10E-6 10E-3])
%     if i == 4
%         cb = colorbar;
%         set(cb,'location','southoutside')
%     end
% end
% 
% %Next plot depth profiles of Uw (wave orbital velocity)
% cc2 = [0 57 255;0 219 218;243 245 4;219 167 0;255 3 3]./255;
% markers = {'o';'s';'d';'^';'p'};
% lineStyle = {'-';'--';'-.';':';'-'};
% spIDs = [2 5 8 11];
% for i = 1:length(hModel)
%     sp(i,2) = subplot(4,3,spIDs(i));
%       
%     %Plot elevation of coral canopy as dashed black line
%     plot(linspace(10E-6,10E-1,10),ones(1,10)*zPatch(i),'--k','linewidth',1.5)
%     hold on
%     
%     %Organize plots by water depth!
%     modelID = find(LR.h == hModel(i) & LR.Hs == 0.4);
%     
%     %Calculate depth profiles
%     for j = 1:length(modelID)
%         Uw = squeeze(LR.Uw(:,:,modelID(j)));
%         xBins = squeeze(LR.xBins(:,:,modelID(j)));
%         zBins = squeeze(LR.zBins(:,:,modelID(j)));
%         [~,xLo] = min(abs(xBins-patchLo));
%         [~,xHi] = min(abs(xBins-patchHi));
%         Profile = nanmean(Uw(xLo:xHi,:),1);
%         z = zBins;
%         x = Profile;
%         plot(x,z,lineStyle{j},...
%             'color',cc2(j,:),'linewidth',1);
%         plot(x(1:7:end),z(1:7:end),markers{j},...
%             'color',cc2(j,:),'markerfacecolor','w',...
%             'markersize',3,'linewidth',1)
%     end
% end
% 
% %Lastly plot depth profiles of normalized turbulence (eps/Uw^3)
% spIDs = [3 6 9 12];
% for i = 1:length(hModel)
%     sp(i,3) = subplot(4,3,spIDs(i));
%       
%     %Plot elevation of coral canopy as dashed black line
%     plot(linspace(10E-6,10E-1,10),ones(1,10)*zPatch(i),'--k','linewidth',1.5)
%     hold on
%     
%     %Organize plots by water depth!
%     modelID = find(LR.h == hModel(i) & LR.Hs == 0.4);
%     
%     %Calculate depth profiles
%     for j = 1:length(modelID)
%         eps = squeeze(LR.eps(:,:,modelID(j)));
%         Uw = squeeze(LR.Uw(:,:,modelID(j)));
%         xBins = squeeze(LR.xBins(:,:,modelID(j)));
%         zBins = squeeze(LR.zBins(:,:,modelID(j)));
%         [~,xLo] = min(abs(xBins-patchLo));
%         [~,xHi] = min(abs(xBins-patchHi));
%         Profile = nanmean(eps(xLo:xHi,:))./nanmean(Uw(xLo:xHi,:).^3);
%         z = zBins;
%         x = Profile;
%         plot(x,z,lineStyle{j},...
%             'color',cc2(j,:),'linewidth',1);
%         plot(x(1:7:end),z(1:7:end),markers{j},...
%             'color',cc2(j,:),'markerfacecolor','w',...
%             'markersize',3,'linewidth',1)
%     end
% end
% wavePeriods = unique(LR.Tp);
% %Make some dummy data for the legends
% xs = repmat(linspace(1E5,1E6,2),length(wavePeriods),1);
% ys = repmat(linspace(min(x),max(x),2),length(wavePeriods),1);
% pp = zeros(length(wavePeriods),1);
% legText = cell(length(wavePeriods),1);
% for i = 1:length(wavePeriods)
%     pp(i) = plot(xs(i,:),ys(i,:),sprintf('%s%s',lineStyle{i},markers{i}),...
%         'markerfacecolor','w','color',cc2(i,:),'linewidth',1);hold on;
%     legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
% end
% leg = legend(pp,legText,'location','southeast','box','on','interpreter','latex');
% 
% %Global adjustments
% set(sp(:,3),'xscale','log','xlim',[10E-6 5E-1],'yticklabel',[])
% set(sp(:,2),'xlim',[0 0.6],'yticklabel',[])
% set(sp(1,:),'ylim',[-1.75 -0.25],'ytick',-1.5:0.5:0,'xticklabel',[])
% set(sp(2,:),'ylim',[-2.75 -1.25],'ytick',-2.5:0.5:-1,'xticklabel',[])
% set(sp(3,:),'ylim',[-3.75 -2.25],'ytick',-3.5:0.5:-2,'xticklabel',[])
% set(sp(4,:),'ylim',[-4.75 -3.25],'ytick',-4.5:0.5:-3)
%  
% %Positoning
% set(sp(1,1),'position',[0.1 0.78 0.25 0.2])
% set(sp(1,2),'position',[0.37 0.78 0.25 0.2])
% set(sp(1,3),'position',[0.64 0.78 0.25 0.2])
% 
% set(sp(2,1),'position',[0.1 0.56 0.25 0.2])
% set(sp(2,2),'position',[0.37 0.56 0.25 0.2])
% set(sp(2,3),'position',[0.64 0.56 0.25 0.2])
% 
% set(sp(3,1),'position',[0.1 0.34 0.25 0.2])
% set(sp(3,2),'position',[0.37 0.34 0.25 0.2])
% set(sp(3,3),'position',[0.64 0.34 0.25 0.2])
% 
% set(sp(4,1),'position',[0.1 0.12 0.25 0.2])
% set(sp(4,2),'position',[0.37 0.12 0.25 0.2])
% set(sp(4,3),'position',[0.64 0.12 0.25 0.2])
% 
% 
% set(leg,'position',[0.86 0.84 0.1 0.1])
% set(cb,'position',[0.1 0.05 0.25 0.015],'linewidth',1)
% 
% %Labeling
% xlabel(cb,'$\overline{\varepsilon} \  \mathrm{(m^2/s^3)}$','interpreter','latex')
% xlabel(sp(4,1),'$\mathrm{Across \ Shore \ Distance \ (m)}$','interpreter','latex')
% xlabel(sp(4,2),'$\langle \overline{U_w}\rangle \  \mathrm{(m/s)}$','interpreter','latex')
% xlabel(sp(4,3),'$\langle \overline{\varepsilon}\rangle/\langle{U_w}^3 \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
% ylabel(sp(1,1),'$z$ (m)','interpreter','latex')
% ylabel(sp(2,1),'$z$ (m)','interpreter','latex')
% ylabel(sp(3,1),'$z$ (m)','interpreter','latex')
% ylabel(sp(4,1),'$z$ (m)','interpreter','latex')
% 
% prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.008 0.008])
%  
% set(ff,'units','inches','renderer','painters');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% % print(ff,[workingDir modelName(1:end-4) '_LRdepthProfiles'],'-dpdf','-r0')
% 
