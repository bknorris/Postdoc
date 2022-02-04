% Create a plot of normalized TKE dissipation for all wave periods at a set
% wave height/water depth. TKE dissipation is normalized by the wave
% orbital velocity (Uw) calculated at every model grid-cell
%
%
% Updates:
% 12/1/21 - Curt asked to compute difference between Low and High Relief.
% Added third panel
%
% This is version 2 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%Define data paths
workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\Figures\';
HRmodelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\HR_Domain\';
LRmodelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\LR_Domain\';
modelInfo = 'Models_allHR_LR.csv';

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Load processed model data
HR = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V4.mat']);
LR = load([LRmodelDir 'LR_postProcess_epsilon_allHR_LR_V4.mat']);

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);
wavePeriods = unique(modelInfo.Tp);
lineStyle = {'-';'--';'-.';':';'-'};
markers = {'o';'s';'d';'^';'p'};
sp = zeros(2,1);
cc = hex2rgb({'#1e5cb3';'#04a1e6';'#ffdf28';'#f58b35';'#bf171f'}); %blue-yellow-red colormap based on GMT panopoly

%Do the LR models first; make sure these are the same as
%Paper1_plotDepthProfiles.m!
patchLo = 40.5;
patchHi = 44;
zPatch = [-0.92 -1.92 -2.92 -3.92];
for i = 1:length(wavePeriods)
    %Organize plots by water depth!
    modelID = find(LR.Tp == wavePeriods(i));
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
    
    sp(1) = subplot(1,3,1);hold on;
    x = LR.gamma(modelID);
    y = epsNorm;
    LR.epsNorm(1:length(modelID),i) = y;
    plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
    B = polyfit(x,log10(y),1);
    Yfit = 10.^(B(1)*linspace(0.1,0.7,length(x))+B(2));
    plot(linspace(0.1,0.7,length(x)),Yfit,'-.',...
        'color',cc(i,:),'linewidth',1)
end

%Then do HR models
patchLo = 35;
patchHi = 41;
zPatch = [-0.9 -1.9 -2.9 -3.9];

for i = 1:length(wavePeriods)
    %Organize plots by water depth!
    modelID = find(HR.Tp == wavePeriods(i));
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
    
    sp(2) = subplot(1,3,2);hold on;
    x = HR.gamma(modelID);
    y = epsNorm;
    HR.epsNorm(1:length(modelID),i) = y;
    plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
    B = polyfit(x,log10(y),1);
    Yfit = 10.^(B(1)*linspace(0.1,0.7,length(x))+B(2));
    plot(linspace(0.1,0.7,length(x)),Yfit,'-.',...
        'color',cc(i,:),'linewidth',1)
end

%Now plot the difference
for i = 1:length(wavePeriods)
    sp(3) = subplot(1,3,3);hold on;
    modelID1 = find(LR.Tp == wavePeriods(i));
    modelID2 = find(HR.Tp == wavePeriods(i));
    if i == 1
        plot([0 25],[0 0],'k','linewidth',1)
    end
    x = LR.gamma(modelID1);
    y = LR.epsNorm(1:length(modelID1),i)-HR.epsNorm(1:length(modelID2),i);
    plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
%     symlog('y')
end

set([sp(1) sp(2)],'yscale','log',...
    'ylim',[5E-3 10E0],...
    'xlim',[0.05 0.75],...
    'xtick',0.1:0.1:0.7)

set(sp(3),'ylim',[-0.5 0.5],...
    'ytick',-0.5:0.25:0.5,...
    'xlim',[0.05 0.75],...
    'xtick',0.1:0.1:0.7)

%Make some dummy data for the legends
xs = repmat(linspace(1E5,1E6,2),length(wavePeriods),1);
ys = repmat(linspace(min(x),max(x),2),length(wavePeriods),1);
pp3 = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    pp3(i) = plot(xs(i,:),ys(i,:),sprintf('%s%s','-.',markers{i}),...
        'markerfacecolor',cc(i,:),'markersize',4,'color',cc(i,:),'linewidth',1);hold on;
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end
leg = legend(pp3,legText,'location','southeast','box','off','interpreter','latex');
 
%Positoning
set(sp(1),'position',[0.1 0.12 0.25 0.8])
set(sp(2),'position',[0.37 0.12 0.25 0.8],'yticklabel',[])
set(sp(3),'position',[0.73 0.12 0.25 0.8])
set(leg,'position',[0.5 0.18 0.1 0.1])

%Labeling
xlabel(sp(1),['$H_s/h \ $' '(-)'],'interpreter','latex')
xlabel(sp(2),['$H_s/h \ $' '(-)'],'interpreter','latex')
xlabel(sp(3),['$H_s/h \ $' '(-)'],'interpreter','latex')
ylabel(sp(1),'$\langle \overline{\varepsilon}\rangle/\langle \overline{{U_w}^3} \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
ylabel(sp(3),['${\langle \overline{\varepsilon}\rangle/\langle \overline{{U_w}^3} \rangle}_{LR}-',...
    '{\langle \overline{\varepsilon}\rangle/\langle \overline{{U_w}^3} \rangle}_{HR}',...
    '\  \mathrm{(m^{-1})}$'],'interpreter','latex')
title(sp(1),'Low Relief')
title(sp(2),'High Relief')
title(sp(3),'Difference')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.005 0.005])

set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 
print(ff,[workingDir 'AllModels_NormalizedDissipation_V2_alt2'],'-dpdf','-r0')