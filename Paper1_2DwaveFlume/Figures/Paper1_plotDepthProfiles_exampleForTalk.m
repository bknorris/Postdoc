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
saveFig = 1;

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
    'position',[100 80   950   400]);

disp('Plotting HR models') %Loop through model scenarios and plot the data
hModel = unique(HR.h);
patchLo = 37;
patchHi = 41.2;
zPatch = [-0.9 -1.9 -2.9 -3.9];
sp = zeros(1,3);
% cc1 = flipud(csvread('c:\Users\user\Documents\Figures\ColorMaps\cequal.csv')./255);
cc1 = hot(265);

sp(1) = subplot(1,3,1);

%Organize plots by water depth!
modelID = find(HR.h == 2 & HR.Hs == 0.4 & HR.Tp == 12);

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
bathyID = find(hModel(2) == hBathy);
data = stlread([HRbathyDir HRmodelBathy{bathyID}]);
[bathyX,sortID] = sort(data.Points(:,1));
bathyY = data.Points(sortID,3);
plot(bathyX-30,bathyY-0.01,'-k','linewidth',1)
set(gca,'ydir','normal','xlim',[patchLo patchHi],'xtick',patchLo:1:patchHi)
colormap(cc1)
caxis([10E-6 10E-3])
cb = colorbar;
set(cb,'location','southoutside')

%Next plot depth profiles of Uw (wave orbital velocity)
cc2 = hex2rgb({'#1e5cb3';'#04a1e6';'#ffdf28';'#f58b35';'#bf171f'}); %blue-yellow-red colormap based on GMT panopoly
% cc2 = hex2rgb({'#1e5cb3';'#1180cd';'#0ea5de';'#7ccb89';'#cce64b'});
% cc2 = flipud(hex2rgb({'#010101';'#0627a3';'#7462fe';'#f58b35';'#fece01'})); %black-purple-orange-yellow colormap based on cequal
% markers = {'^';'s';'d';'o';'p'};
lineStyle = {'-';'--';'-.';':';'-'};
sp(2) = subplot(1,3,2);

%Plot elevation of coral canopy as dashed black line
plot(linspace(10E-6,10E-1,10),ones(1,10)*zPatch(2),'--k','linewidth',1.5)
hold on

%Organize plots by water depth!
modelID = find(HR.h == 2 & HR.Hs == 0.4);

%Calculate depth profiles
for j = 1:length(modelID)
    Uw = squeeze(HR.Uw(:,:,modelID(j)));
    xBins = squeeze(HR.xBins(:,:,modelID(j)));
    zBins = squeeze(HR.zBins(:,:,modelID(j)));
    [~,xLo] = min(abs(xBins-patchLo));
    [~,xHi] = min(abs(xBins-patchHi));
    Profile = nanmean(Uw(xLo:xHi,:),1);
    Profile(1:5) = NaN;
    z = zBins;
    x = Profile;
    plot(x,z,lineStyle{j},...
        'color',cc2(j,:),'linewidth',1.5);
    %         plot(x(1:7:end),z(1:7:end),markers{j},...
    %             'color',cc2(j,:),'markerfacecolor',cc2(j,:),...
    %             'markersize',3,'linewidth',1)
end

%Lastly plot depth profiles of normalized turbulence (eps/Uw^3)
sp(3) = subplot(1,3,3);

%Plot elevation of coral canopy as dashed black line
plot(linspace(10E-6,10E-1,10),ones(1,10)*zPatch(2),'--k','linewidth',1.5)
hold on

%Organize plots by water depth!
modelID = find(HR.h == 2 & HR.Hs == 0.4);

%Calculate depth profiles
for j = 1:length(modelID)
    eps = squeeze(HR.eps(:,:,modelID(j)));
    Uw = squeeze(HR.Uw(:,:,modelID(j)));
    xBins = squeeze(HR.xBins(:,:,modelID(j)));
    zBins = squeeze(HR.zBins(:,:,modelID(j)));
    [~,xLo] = min(abs(xBins-patchLo));
    [~,xHi] = min(abs(xBins-patchHi));
    Profile = nanmean(eps(xLo:xHi,:))./nanmean(Uw(xLo:xHi,:).^3);
    Profile(1:5) = NaN;
    z = zBins;
    x = Profile;
    plot(x,z,lineStyle{j},...
        'color',cc2(j,:),'linewidth',1.5);
    %         plot(x(1:7:end),z(1:7:end),markers{j},...
    %             'color',cc2(j,:),'markerfacecolor',cc2(j,:),...
    %             'markersize',3,'linewidth',1)
end

wavePeriods = unique(HR.Tp);
%Make some dummy data for the legends
xs = repmat(linspace(1E5,1E6,2),length(wavePeriods),1);
ys = repmat(linspace(min(x),max(x),2),length(wavePeriods),1);
pp = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    pp(i) = plot(xs(i,:),ys(i,:),sprintf('%s',lineStyle{i}),...
        'color',cc2(i,:),'linewidth',1.5);hold on;
%     pp(i) = plot(xs(i,:),ys(i,:),sprintf('%s%s',lineStyle{i},markers{i}),...
%         'markerfacecolor',cc2(i,:),'color',cc2(i,:),'linewidth',1);hold on;
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end
leg = legend(pp,legText,'location','southeast','box','on','interpreter','latex');

%Global adjustments
set(sp(3),'xscale','log','xlim',[10E-6 5E-1],'yticklabel',[])
set(sp(2),'xlim',[0 0.6],'yticklabel',[])
set(sp,'ylim',[-2.75 -1.25],'ytick',-2.5:0.5:-1)

%Positoning
set(sp(1),'position',[0.1 0.25 0.25 0.72])
set(sp(2),'position',[0.37 0.25 0.25 0.72])
set(sp(3),'position',[0.64 0.25 0.25 0.72])

set(leg,'position',[0.84 0.75 0.1 0.1])
set(cb,'position',[0.1 0.1 0.25 0.025],'linewidth',1)

%Labeling
xlabel(cb,'$\overline{\varepsilon} \  \mathrm{(m^2/s^3)}$','interpreter','latex')
xlabel(sp(1),'$\mathrm{Across \ Shore \ Distance \ (m)}$','interpreter','latex')
xlabel(sp(2),'$\langle \overline{U_w}\rangle \  \mathrm{(m/s)}$','interpreter','latex')
xlabel(sp(3),'$\langle \overline{\varepsilon}\rangle/\langle \overline{{U_w}^3} \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
ylabel(sp(1),'$z$ (m)','interpreter','latex')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.008 0.008])

set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

if saveFig == 1
    print(ff,[workingDir modelName(1:end-4) '_HRdepthProfileForTalk'],'-dpdf','-r0')
end
%% Create Figure
ff = figure(2);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   950   400]);

disp('Plotting LR models') %Loop through model scenarios and plot the data
hModel = unique(LR.h);
patchLo = 39;
patchHi = 43.2;
zPatch = [-1.25 -2.25 -3.25 -4.25];
sp = zeros(1,3);
% cc1 = flipud(csvread('c:\Users\user\Documents\Figures\ColorMaps\cequal.csv')./255);
cc1 = hot(265);

sp(1) = subplot(1,3,1);

%Organize plots by water depth!
modelID = find(LR.h == 2 & LR.Hs == 0.4 & LR.Tp == 12);

xBins = squeeze(LR.xBins(:,:,modelID));
zBins = squeeze(LR.zBins(:,:,modelID));
[~,xLo] = min(abs(xBins-patchLo));
[~,xHi] = min(abs(xBins-patchHi));
eps = squeeze(LR.eps(:,:,modelID));
pcolor(xBins(xLo:xHi),zBins,eps(xLo:xHi,:)')
shading interp
set(gca,'ydir','normal')
hold on
%Now plot the bathy as a black line
hBathy = regexp(LRmodelBathy,'(?<=_)(.*?)(?=m)','match');
hBathy = str2double([hBathy{:}]);
bathyID = find(hModel(2) == hBathy);
data = stlread([LRbathyDir LRmodelBathy{bathyID}]);
[bathyX,sortID] = sort(data.Points(:,1));
bathyY = data.Points(sortID,3);
plot(bathyX-30,bathyY-0.01,'-k','linewidth',1)
set(gca,'ydir','normal','xlim',[patchLo patchHi],'xtick',patchLo:1:patchHi)
colormap(cc1)
caxis([10E-6 10E-3])
cb = colorbar;
set(cb,'location','southoutside')

%Next plot depth profiles of Uw (wave orbital velocity)
cc2 = hex2rgb({'#1e5cb3';'#04a1e6';'#ffdf28';'#f58b35';'#bf171f'}); %blue-yellow-red colormap based on GMT panopoly
% cc2 = hex2rgb({'#1e5cb3';'#1180cd';'#0ea5de';'#7ccb89';'#cce64b'});
% cc2 = flipud(hex2rgb({'#010101';'#0627a3';'#7462fe';'#f58b35';'#fece01'})); %black-purple-orange-yellow colormap based on cequal
% markers = {'^';'s';'d';'o';'p'};
lineStyle = {'-';'--';'-.';':';'-'};
sp(2) = subplot(1,3,2);

%Plot elevation of coral canopy as dashed black line
plot(linspace(10E-6,10E-1,10),ones(1,10)*zPatch(2),'--k','linewidth',1.5)
hold on

%Organize plots by water depth!
modelID = find(LR.h == 2 & LR.Hs == 0.4);

%Calculate depth profiles
for j = 1:length(modelID)
    Uw = squeeze(LR.Uw(:,:,modelID(j)));
    xBins = squeeze(LR.xBins(:,:,modelID(j)));
    zBins = squeeze(LR.zBins(:,:,modelID(j)));
    [~,xLo] = min(abs(xBins-patchLo));
    [~,xHi] = min(abs(xBins-patchHi));
    Profile = nanmean(Uw(xLo:xHi,:),1);
    Profile(1:5) = NaN;
    z = zBins;
    x = Profile;
    plot(x,z,lineStyle{j},...
        'color',cc2(j,:),'linewidth',1.5);
    %         plot(x(1:7:end),z(1:7:end),markers{j},...
    %             'color',cc2(j,:),'markerfacecolor',cc2(j,:),...
    %             'markersize',3,'linewidth',1)
end

%Lastly plot depth profiles of normalized turbulence (eps/Uw^3)
sp(3) = subplot(1,3,3);

%Plot elevation of coral canopy as dashed black line
plot(linspace(10E-6,10E-1,10),ones(1,10)*zPatch(2),'--k','linewidth',1.5)
hold on

%Organize plots by water depth!
modelID = find(LR.h == 2 & LR.Hs == 0.4);

%Calculate depth profiles
for j = 1:length(modelID)
    eps = squeeze(LR.eps(:,:,modelID(j)));
    Uw = squeeze(LR.Uw(:,:,modelID(j)));
    xBins = squeeze(LR.xBins(:,:,modelID(j)));
    zBins = squeeze(LR.zBins(:,:,modelID(j)));
    [~,xLo] = min(abs(xBins-patchLo));
    [~,xHi] = min(abs(xBins-patchHi));
    Profile = nanmean(eps(xLo:xHi,:))./nanmean(Uw(xLo:xHi,:).^3);
    Profile(1:5) = NaN;
    z = zBins;
    x = Profile;
    plot(x,z,lineStyle{j},...
        'color',cc2(j,:),'linewidth',1.5);
    %         plot(x(1:7:end),z(1:7:end),markers{j},...
    %             'color',cc2(j,:),'markerfacecolor',cc2(j,:),...
    %             'markersize',3,'linewidth',1)
end

wavePeriods = unique(LR.Tp);
%Make some dummy data for the legends
xs = repmat(linspace(1E5,1E6,2),length(wavePeriods),1);
ys = repmat(linspace(min(x),max(x),2),length(wavePeriods),1);
pp = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    pp(i) = plot(xs(i,:),ys(i,:),sprintf('%s',lineStyle{i}),...
        'color',cc2(i,:),'linewidth',1.5);hold on;
%     pp(i) = plot(xs(i,:),ys(i,:),sprintf('%s%s',lineStyle{i},markers{i}),...
%         'markerfacecolor',cc2(i,:),'color',cc2(i,:),'linewidth',1);hold on;
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end
leg = legend(pp,legText,'location','southeast','box','on','interpreter','latex');

%Global adjustments
set(sp(3),'xscale','log','xlim',[10E-6 5E-1],'yticklabel',[])
set(sp(2),'xlim',[0 0.6],'yticklabel',[])
set(sp,'ylim',[-2.75 -1.25],'ytick',-2.5:0.5:-1)

%Positoning
set(sp(1),'position',[0.1 0.25 0.25 0.72])
set(sp(2),'position',[0.37 0.25 0.25 0.72])
set(sp(3),'position',[0.64 0.25 0.25 0.72])

set(leg,'position',[0.84 0.75 0.1 0.1])
set(cb,'position',[0.1 0.1 0.25 0.025],'linewidth',1)

%Labeling
xlabel(cb,'$\overline{\varepsilon} \  \mathrm{(m^2/s^3)}$','interpreter','latex')
xlabel(sp(1),'$\mathrm{Across \ Shore \ Distance \ (m)}$','interpreter','latex')
xlabel(sp(2),'$\langle \overline{U_w}\rangle \  \mathrm{(m/s)}$','interpreter','latex')
xlabel(sp(3),'$\langle \overline{\varepsilon}\rangle/\langle \overline{{U_w}^3} \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
ylabel(sp(1),'$z$ (m)','interpreter','latex')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.008 0.008])

set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

if saveFig == 1
    print(ff,[workingDir modelName(1:end-4) '_LRdepthProfileForTalk'],'-dpdf','-r0')
end

