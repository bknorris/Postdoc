% Create a plot of wave energy flux spectra for an example 4s and 12s model
% from the HR and LR domains, respectively. 
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

modelInfo = 'Models_waveEnergyFlux_example.csv';
% modelDepthScenarios = {'0.5 m';'1 m';'2 m';'3 m';'4 m'};

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Get modelNames from the model directory
HRmodelFree = dir([HRmodelDir '*freeSurf*']);
HRmodelFree = HRmodelFree(~ismember({HRmodelFree.name},{'.','..'}));
HRmodelFree = {HRmodelFree.name};

LRmodelFree = dir([LRmodelDir '*freeSurf_V1*']);
LRmodelFree = LRmodelFree(~ismember({LRmodelFree.name},{'.','..'}));
LRmodelFree = {LRmodelFree.name};

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);
pp1 = length(modelInfo.scenarioNumber);
pp2 = length(modelInfo.scenarioNumber);
cc = [0.6471 0 0.1490;...
    0.1922 0.2118 0.5843];
disp('Loading HR models')
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(HRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),HRmodelFree{modelID})
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
    F1 = 1025*9.81*WG1(2:end).*Cg(2:end);
    F2 = 1025*9.81*WG2(2:end).*Cg(2:end);
    sp(1) = subplot(121);
    pp1(i) = plot(F(2:end),F1./1000,'-','linewidth',1.5,'color',cc(i,:));hold on %convert to KW
    pp2(i) = plot(F(2:end),F2./1000,'--','linewidth',1.5,'color',cc(i,:));
end

pp3 = length(modelInfo.scenarioNumber);
pp4 = length(modelInfo.scenarioNumber);
disp('Loading LR models')
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(LRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),LRmodelFree{modelID})
    load([LRmodelDir LRmodelFree{modelID}]);
    FS = data;clear data
    
    %Calculate wave energy dissipation @ WG1 and WG2 
    [WG1,~] = pwelch(FS.WG1.eta,[],[],[],8);
    [WG2,F] = pwelch(FS.WG2.eta,[],[],[],8);
    
    k=wavek(F,modelInfo.h(i));
    kh = k*modelInfo.h(i);
    
    c = sqrt((9.81/k)*tanh(kh));
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c; 
    F1 = 1025*9.81*WG1(2:end).*Cg(2:end);
    F2 = 1025*9.81*WG2(2:end).*Cg(2:end);
    sp(2) = subplot(122);
    pp3(i) = plot(F(2:end),F1./1000,'-','linewidth',1.5,'color',cc(i,:));hold on
    pp4(i) = plot(F(2:end),F2./1000,'--','linewidth',1.5,'color',cc(i,:));
end
%Global Adjustments
set(sp,'yscale','log','xlim',[0 1])

% leg1 = legend(sp(1),[pp1(1) pp2(1)],{'WG1';'WG2'});
% leg2 = legend(sp(2),[pp3(1) pp4(1)],{'WG1';'WG2'});

set(sp(1),'position',[0.1 0.15 0.4 0.8])
set(sp(2),'position',[0.55 0.15 0.4 0.8],'yticklabel',[])

%Labeling
xlabel(sp(1),['$f \quad$' '(Hz)'],'interpreter','latex')
xlabel(sp(2),['$f \quad$' '(Hz)'],'interpreter','latex')
ylabel(sp(1),'$F_0 \quad (KW/m)$','interpreter','latex')
title(sp(1),'HR Models')
title(sp(2),'LR Models')

prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.01 0.01])

%Save figure
export_fig(ff,[workingDir  'HR_LR_Example_WaveEnergyFluxes'],'-png','-r600','-nocrop')