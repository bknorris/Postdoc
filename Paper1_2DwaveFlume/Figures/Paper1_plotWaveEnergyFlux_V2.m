% Create a plot of spectral wave energy flux for all wave
% periods at a set wave height
%
% Updates:
% 08/19/21 - Modified script to plot wave energy flux for
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
modelInfo = 'Models_waveEnergyFlux.csv';
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

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);
cc = [0 57 255;54 94 228;202 147 147;250 114 114;255 3 3]./255;
disp('Loading HR models')
%Loop through model scenarios and plot the data
pp1 = zeros(length(modelInfo.scenarioNumber),1);
pp2 = zeros(length(modelInfo.scenarioNumber),1);
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(HRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),HRmodelFree{modelID})
    FS = load([HRmodelDir HRmodelFree{modelID}]);
    eta1 = repmat(FS.WG1.eta,1,5);
    eta2 = repmat(FS.WG7.eta,1,5);
    %Calculate wave energy dissipation @ WG1 and WG2 
    [WG1,f] = pwelch(eta1,round(length(eta1)/20,0),[],[],8);
    [WG2,~] = pwelch(eta2,round(length(eta2)/20,0),[],[],8);
    
    k=wavek(f,modelInfo.h(i));
    kh = k*modelInfo.h(i);
    c = sqrt(9.81.*tanh(kh))./sqrt(k);
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c; 
    F1 = 1025*9.81*WG1(2:end).*Cg(2:end);
    F2 = 1025*9.81*WG2(2:end).*Cg(2:end);
    sp(1) = subplot(121);
    pp1(i) = plot(f(2:end),F1,'-','linewidth',1.5,'color',cc(i,:));hold on;
    pp2(i) = plot(f(2:end),F2,'--','linewidth',1.5,'color',cc(i,:));
end
disp('Loading LR models')
%Loop through model scenarios and plot the data
pp3 = length(modelInfo.scenarioNumber);
pp4 = length(modelInfo.scenarioNumber);
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(LRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),LRmodelFree{modelID})
    FS = load([LRmodelDir LRmodelFree{modelID}]);
    eta1 = repmat(FS.WG1.eta,1,5);
    eta2 = repmat(FS.WG7.eta,1,5);
    %Calculate wave energy dissipation @ WG1 and WG2 
    [WG1,f] = pwelch(eta1,round(length(eta1)/20,0),[],[],8);
    [WG2,~] = pwelch(eta2,round(length(eta2)/20,0),[],[],8);
    
    k=wavek(f,modelInfo.h(i));
    kh = k*modelInfo.h(i);
    c = sqrt(9.81.*tanh(kh))./sqrt(k);
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c; 
    F1 = 1025*9.81*WG1(2:end).*Cg(2:end);
    F2 = 1025*9.81*WG2(2:end).*Cg(2:end);
    sp(2) = subplot(122);
    pp3(i) = plot(f(2:end),F1,'-','linewidth',1.5,'color',cc(i,:));hold on;
    pp4(i) = plot(f(2:end),F2,'--','linewidth',1.5,'color',cc(i,:));
end
%Global Adjustments
set(sp,'yscale','log','xlim',[0 1])

set(sp(1),'position',[0.1 0.15 0.4 0.8])
set(sp(2),'position',[0.55 0.15 0.4 0.8],'yticklabel',[])

%Labeling
xlabel(sp(1),['$f \ $' '(Hz)'],'interpreter','latex')
xlabel(sp(2),['$f \ $' '(Hz)'],'interpreter','latex')
ylabel(sp(1),'$F_0 \ $ (W/m)','interpreter','latex')
title(sp(1),'HR Models')
title(sp(2),'LR Models')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.005 0.005])
 
set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% print(ff,[workingDir 'AllModels_WaveEnergyDissipation'],'-dpdf','-r0')