% Create a plot of normalized TKE dissipation for all wave periods at a set
% wave height/water depth. TKE dissipation is normalized by the wave
% orbital velocity (Uw) calculated at every model grid-cell
%
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
HR = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V5.mat']);
LR = load([LRmodelDir 'LR_postProcess_epsilon_allHR_LR_V5.mat']);

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);
wavePeriods = unique(modelInfo.Tp);
waterDepth = unique(modelInfo.h);
lineStyle = {'-';'--';'-.';':';'-'};
markers = {'o';'s';'d';'^';'p'};
sp = zeros(2,1);
cc = hex2rgb({'#1e5cb3';'#04a1e6';'#ffdf28';'#f58b35';'#bf171f'}); %blue-yellow-red colormap based on GMT panopoly
% cc = flipud(hex2rgb({'#032760';'#0057d2';'#57b0ff';'#73dcff'})); %A dark -> light blue colormap (based on Bath 112)
%Do the HR models first; make sure these are the same as
%Paper1_plotDepthProfiles.m!
patchLo = 35;
patchHi = 41;
zPatch = [-0.9 -1.9 -2.9 -3.9];
upper = 15;
for i = 1:length(wavePeriods)
    %Collect data by wave period then water depth!
    modelID = find(HR.Tp == wavePeriods(i));
    eps = zeros(length(modelID),1);
    Uw = zeros(length(modelID),1);
    
    %EXPERIMENTAL: try filtering results by their proximity to the bed
    for j = 1:length(modelID)
        eps1 = squeeze(HR.eps(:,:,modelID(j)));
        Uw1 = squeeze(HR.Ubr(:,:,modelID(j)));
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
        
        %Find the "bed" as the last value in the column that is a nan
        epscrop = eps1(xLo:xHi,:);[m,n] = size(epscrop);
        uwcrop = Uw1(xLo:xHi,:);
        epsfilt = NaN(m,n);
        uwfilt = NaN(m,n);
        for k = 1:m
            nanID = find(isnan(epscrop(k,1:end-1)),1,'last');
            epsfilt(k,nanID+1:nanID+upper) = epscrop(k,nanID+1:nanID+upper);
            nanID = find(isnan(uwcrop(k,1:end-1)),1,'last');
            uwfilt(k,nanID+1:nanID+upper) = uwcrop(k,nanID+1:nanID+upper);
        end
            
        eps(j) = nanmean(nanmean(epsfilt));
        Uw(j) = nanmean(nanmean(uwfilt.^3));
    end
    sp(1) = subplot(1,2,1);hold on;
    x = Uw;
    y = eps;
    plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
    %calculate wave friction factor (fw)
    B = polyfit(x,y,1);
    fw = B(1)./(2/(3*pi));
%     fw = mean((3*pi.*eps)./(2*Uw));
    fprintf('Tp = %0.2f\n',wavePeriods(i))
    fprintf('fw = %0.3f\n\n',fw)
end

%Then do LR models
patchLo = 40.5;
patchHi = 43.5;
zPatch = [-0.92 -1.92 -2.92 -3.92];

for i = 1:length(wavePeriods)
    %Collect data by wave period then water depth!
    modelID = find(LR.Tp == wavePeriods(i));
    eps = zeros(length(modelID),1);
    Uw = zeros(length(modelID),1);
    
    %EXPERIMENTAL: try filtering results by their proximity to the bed
    for j = 1:length(modelID)
        eps1 = squeeze(LR.eps(:,:,modelID(j)));
        Uw1 = squeeze(LR.Ubr(:,:,modelID(j)));
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
        
        %Find the "bed" as the last value in the column that is a nan
        epscrop = eps1(xLo:xHi,:);[m,n] = size(epscrop);
        uwcrop = Uw1(xLo:xHi,:);
        epsfilt = NaN(m,n);
        uwfilt = NaN(m,n);
        for k = 1:m
            nanID = find(isnan(epscrop(k,1:end-1)),1,'last');
            epsfilt(k,nanID+1:nanID+upper) = epscrop(k,nanID+1:nanID+upper);
            nanID = find(isnan(uwcrop(k,1:end-1)),1,'last');
            uwfilt(k,nanID+1:nanID+upper) = uwcrop(k,nanID+1:nanID+upper);
        end
            
        eps(j) = nanmean(nanmean(epsfilt));
        Uw(j) = nanmean(nanmean(uwfilt.^3));
    end
    sp(2) = subplot(1,2,2);hold on;
    x = Uw;
    y = eps;
    plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
    %calculate wave friction factor (fw)
    B = polyfit(x,y,1);
    fw = B(1)./(2/(3*pi));
%     fw = mean((3*pi.*eps)./(2*Uw));
    fprintf('Tp = %0.2f\n',wavePeriods(i))
    fprintf('fw = %0.3f\n\n',fw)
end

%Make some dummy data for the legends
xs = repmat(linspace(1E5,1E6,2),length(wavePeriods),1);
ys = repmat(linspace(min(x),max(x),2),length(wavePeriods),1);
pp3 = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    pp3(i) = plot(xs(i,:),ys(i,:),sprintf('%s',markers{i}),...
        'markerfacecolor',cc(i,:),'markersize',4,'color',cc(i,:),'linewidth',1);hold on;
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end
leg = legend(pp3,legText,'location','southeast','box','off','interpreter','latex');

%Axes and plot adjustments
set(sp,'yscale','log',...
    'xscale','log',...
    'ylim',[10E-7 10E-2],...
    'xlim',[10E-4 5])
set(sp(1),'position',[0.1 0.15 0.35 0.8])
set(sp(2),'position',[0.5 0.15 0.35 0.8],'yticklabel',[])
set(leg,'position',[0.73 0.18 0.1 0.2])

%Labeling
xlabel(sp(1),['$H_s/h \ $' '(-)'],'interpreter','latex')
xlabel(sp(2),['$H_s/h \ $' '(-)'],'interpreter','latex')
xlabel(sp(1),'$\langle \overline{{U_w}^3} \rangle \  \mathrm{(m^3/s^3)}$','interpreter','latex')
xlabel(sp(2),'$\langle \overline{{U_w}^3} \rangle \  \mathrm{(m^3/s^3)}$','interpreter','latex')
ylabel(sp(1),'$\langle \overline{\varepsilon} \rangle \  \mathrm{(m^2/s^3)}$','interpreter','latex')
title(sp(1),'HR Models')
title(sp(2),'LR Models')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.005 0.005])

set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% print(ff,[workingDir 'AllModels_EpsVsUcubed'],'-dpdf','-r0')