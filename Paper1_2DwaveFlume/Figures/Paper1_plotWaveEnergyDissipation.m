% Create a plot of spectral wave energy dissipation spectra for all wave
% periods at a set wave height
%
% Updates:
% 08/19/21 - Modified script to plot wave energy dissipation (not flux) for
%            five models (Tp = 4, 8, 12, 16, 20 s). 
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
% modelDepthScenarios = {'0.5 m';'1 m';'2 m';'3 m';'4 m'};

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Get modelNames from the model directory
HRmodelFree = dir([HRmodelDir '*freeSurf_V1*']);
HRmodelFree = HRmodelFree(~ismember({HRmodelFree.name},{'.','..'}));
HRmodelFree = {HRmodelFree.name};

LRmodelFree = dir([LRmodelDir '*freeSurf_V1*']);
LRmodelFree = LRmodelFree(~ismember({LRmodelFree.name},{'.','..'}));
LRmodelFree = {LRmodelFree.name};

disp('Loading HR models')
%Loop through model scenarios and plot the data
DHR = zeros(length(modelInfo.scenarioNumber),1);
GammaHR = zeros(length(modelInfo.scenarioNumber),1);
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(HRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),HRmodelFree{modelID})
    FS = load([HRmodelDir HRmodelFree{modelID}]);
    
    %Calculate wave energy dissipation @ WG1 and WG2 
    [WG1,~] = pwelch(FS.WG1.eta,[],[],[],8);
    [WG2,F] = pwelch(FS.WG2.eta,[],[],[],8);
    
    k=wavek(F,modelInfo.h(i));
    kh = k*modelInfo.h(i);
    c = sqrt(9.81.*tanh(kh))./sqrt(k);
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c;
    cutoff = find(F>=0.5,1,'last');
    F1 = trapz(1025*9.81*WG1(2:cutoff).*Cg(2:cutoff));
    F2 = trapz(1025*9.81*WG2(2:cutoff).*Cg(2:cutoff));
    DHR(i) = abs((F1-F2)/(73-31));
    GammaHR(i) = modelInfo.Hs(i)/modelInfo.h(i);
end
disp('Loading LR models')
%Loop through model scenarios and plot the data
DLR = zeros(length(modelInfo.scenarioNumber),1);
GammaLR = zeros(length(modelInfo.scenarioNumber),1);
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(LRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),LRmodelFree{modelID})
    FS = load([LRmodelDir LRmodelFree{modelID}]);
    
    %Calculate wave energy dissipation @ WG1 and WG2 
    [WG1,~] = pwelch(FS.WG1.eta,[],[],[],8);
    [WG2,F] = pwelch(FS.WG2.eta,[],[],[],8);
    
    k=wavek(F,modelInfo.h(i));
    kh = k*modelInfo.h(i);
    
    c = sqrt((9.81/k)*tanh(kh));
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c; 
    cutoff = find(F>=0.5,1,'last');
    F1 = trapz(1025*9.81*WG1(2:cutoff).*Cg(2:cutoff));
    F2 = trapz(1025*9.81*WG2(2:cutoff).*Cg(2:cutoff));
    DLR(i) = abs((F1-F2)/(73-31));
    GammaLR(i) = modelInfo.Hs(i)/modelInfo.h(i);
end

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);
wavePeriods = unique(modelInfo.Tp);
waterDepths = unique(modelInfo.h);
lineStyle = {'-';'--';'-.';':';'-'};
markers = {'o';'s';'d';'^';'p'};
pp1 = zeros(length(wavePeriods),1);
pp2 = zeros(length(wavePeriods),1);
cc = [0 57 255;54 94 228;202 147 147;250 114 114;255 3 3]./255;
for i = 1:length(wavePeriods)
    Tp = wavePeriods(i);
    idx = find(modelInfo.Tp == wavePeriods(i));
    sp(1) = subplot(1,2,1);hold on;
    x = GammaHR(idx);y = DHR(idx);
    pp1(i) = plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor','w','linewidth',1);
    B = polyfit(x,log10(y),1);
    Yfit = 10.^(B(1)*linspace(min(x),max(x),length(x))+B(2));
    plot(linspace(min(x),max(x),length(x)),Yfit,sprintf('%s',lineStyle{i}),...
        'color',cc(i,:),'markerfacecolor','w','linewidth',1)
    sp(2) = subplot(1,2,2);hold on;
    x = GammaLR(idx);y = DLR(idx);
    pp2(i) = plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor','w','linewidth',1);
    B = polyfit(x,log10(y),1);
    Yfit = 10.^(B(1)*linspace(min(x),max(x),length(x))+B(2));
    plot(linspace(min(x),max(x),length(x)),Yfit,sprintf('%s',lineStyle{i}),...
        'color',cc(i,:),'markerfacecolor','w','linewidth',1)
end
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

%Axes and plot adjustments
set(sp,'yscale','log',...
    'ylim',[1 20000],...
    'xlim',[0.05 0.7],...
    'xtick',0.1:0.1:0.7)
set(sp(1),'position',[0.1 0.15 0.35 0.8])
set(sp(2),'position',[0.5 0.15 0.35 0.8],'yticklabel',[])
set(leg,'position',[0.87 0.45 0.1 0.2])

%Labeling
xlabel(sp(1),['$H_s/h \ $' '(-)'],'interpreter','latex')
xlabel(sp(2),['$H_s/h \ $' '(-)'],'interpreter','latex')
ylabel(sp(1),'$\varepsilon_0 \  \mathrm{(W/m^2)}$','interpreter','latex')
title(sp(1),'HR Models')
title(sp(2),'LR Models')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.005 0.005])
 
set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% print(ff,[workingDir 'AllModels_WaveEnergyDissipation'],'-dpdf','-r0')