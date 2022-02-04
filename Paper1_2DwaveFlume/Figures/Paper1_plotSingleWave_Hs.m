% Plot water level time series for a single wave and the significant wave
% height for each of five water depth models specified in a CSV for the HR
% and LR domains. Single waves are plotted from WG4 (in the center of the
% hi-res patches)
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%Define data paths
workingDir = 'd:\USGS\Models\Figures\Paper1_2DwaveFlume\';
HRmodelDir = 'd:\USGS\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\DataAnalysis\';
LRmodelDir = 'd:\USGS\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\DataAnalysis\';
modelInfo = 'Models_Hs0_4m_Tp_4_ramptime8.csv';
modelDepthScenarios = {'0.5 m';'1 m';'2 m';'3 m';'4 m'};
waveGaugeNames = {'WG1';'WG2';'WG3';'WG4';'WG5';'WG6'};
HRwaveGauges = [66 68 70 72 74 76];
LRwaveGauges = [63 65 67 69 71 73];

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Get modelNames from the model directory
HRmodelScenarios = dir([HRmodelDir '*freeSurf*']);
HRmodelScenarios = HRmodelScenarios(~ismember({HRmodelScenarios.name},{'.','..'}));
HRmodelScenarios = {HRmodelScenarios.name};

LRmodelScenarios = dir([LRmodelDir '*freeSurf*']);
LRmodelScenarios = LRmodelScenarios(~ismember({LRmodelScenarios.name},{'.','..'}));
LRmodelScenarios = {LRmodelScenarios.name};

%Create figure
ff = figure(1);
    set(ff,'PaperOrientation','landscape',...
    'position',[100 80   850   500]);
dataToPlot = 66:160; %adjust depending on wave period!
nfft = 45;
pp1 = zeros(5,1);
pp2 = zeros(5,1);
pp3 = zeros(5,1);
pp4 = zeros(5,1);
cc = brewermap(7,'Blues');
disp('Loading HR models')
count = 1;
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(HRmodelScenarios','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    if isempty(modelID) %process only the folders listed in the CSV file
        continue
    else
        fprintf('Loading file %0.0f of %0.0f: %s\n',count,length(modelInfo.scenarioNumber),HRmodelScenarios{modelID})
        load([HRmodelDir HRmodelScenarios{modelID}]);
        HR = data;clear data
        
        %Plot free surface data
        sp(1) = subplot(221);
        pp1(i) = plot(linspace(0,max(dataToPlot)-min(dataToPlot),length(dataToPlot)),HR.WG4.eta(dataToPlot),...
            'color',cc(count+2,:),'linewidth',1.5);
        hold on
        
        sp(2) = subplot(223);
        %Calculate Hs at each WG for each model
        fn = fieldnames(HR);
        Hs = zeros(length(fn),1);
        for j = 1:length(fn)
            P = repmat(HR.(fn{j}).eta(dataToPlot),1,5);
            
            %Spectral analysis
            dt = 1/0.125;
            [Cpp,F] = pwelch(P,[],[],[],dt);
            lfc = find(F>=0.033,1,'first');hfc = find(F<=0.33,1,'last');
            Hs(j) = 4*sqrt(sum(Cpp(lfc:hfc).*mean(diff(F(lfc:hfc))))); %Hs from the model
        end
        plot(linspace(30,50,10),ones(1,10).*modelInfo.Hs,'--k','linewidth',1.5)
        pp2(i) = plot(HRwaveGauges-30,Hs,'o','color',cc(count+2,:),'markerfacecolor',cc(count+2,:),'markersize',6);
        hold on
        count = count+1;
    end
end

disp('Loading LR models')
count = 1;
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(LRmodelScenarios','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    if isempty(modelID) %process only the folders listed in the CSV file
        continue
    else
        fprintf('Loading file %0.0f of %0.0f: %s\n',count,length(modelInfo.scenarioNumber),LRmodelScenarios{modelID})
        load([LRmodelDir LRmodelScenarios{modelID}]);
        LR = data;clear data
        
        %Plot free surface data
        sp(3) = subplot(222);
        pp3(i) = plot(linspace(0,max(dataToPlot)-min(dataToPlot),length(dataToPlot)),LR.WG4.eta(dataToPlot),...
            'color',cc(count+2,:),'linewidth',1.5);
        hold on
        
        sp(4) = subplot(224);
        %Calculate Hs at each WG for each model
        fn = fieldnames(LR);
        Hs = zeros(length(fn),1);
        for j = 1:length(fn)
            P = repmat(LR.(fn{j}).eta(dataToPlot),1,5);
            
            %Spectral analysis
            dt = 1/0.125;
            [Cpp,F] = pwelch(P,[],[],[],dt);
            lfc = find(F>=0.033,1,'first');hfc = find(F<=0.33,1,'last');
            Hs(j) = 4*sqrt(sum(Cpp(lfc:hfc).*mean(diff(F(lfc:hfc))))); %Hs from the model
        end
        plot(linspace(30,50,10),ones(1,10).*modelInfo.Hs,'--k','linewidth',1.5)
        pp4(i) = plot(LRwaveGauges-30,Hs,'o','color',cc(count+2,:),'markerfacecolor',cc(count+2,:),'markersize',6);
        hold on
        count = count+1;
    end
end
prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.005 0.005])

%Plot wave gauges on lower subplots
plot(sp(2),HRwaveGauges-30,repmat(0.65,1,6),'+k','markersize',8);
plot(sp(2),HRwaveGauges-30,repmat(0.65,1,6),'ok','markersize',4);
text(sp(2),HRwaveGauges-30.5,repmat(0.71,1,6),waveGaugeNames,'fontsize',6)

plot(sp(4),LRwaveGauges-30,repmat(0.65,1,6),'+k','markersize',8);
plot(sp(4),LRwaveGauges-30,repmat(0.65,1,6),'ok','markersize',4);
text(sp(4),LRwaveGauges-30.5,repmat(0.71,1,6),waveGaugeNames,'fontsize',6)

%Create legend
leg1 = legend(pp3,modelDepthScenarios);
leg2 = legend(pp4,modelDepthScenarios);
set(leg1,'position',[0.88 0.81 0.05 0.05])
set(leg2,'position',[0.88 0.34 0.05 0.05])

%Global Adjustments
set([sp(1) sp(3)],'xlim',[50 80])
set([sp(2) sp(4)],'ylim',[0 0.8])
set(sp(1),'position',[0.1 0.58 0.35 0.36])
set(sp(2),'position',[0.1 0.1 0.35 0.36])
set(sp(3),'position',[0.5 0.58 0.35 0.36],...
    'yticklabel',[])
set(sp(4),'position',[0.5 0.1 0.35 0.36],...
    'yticklabel',[])

%Labeling
title(sp(1),'HR Models')
title(sp(3),'LR Models')
ylabel(sp(1),'\eta (m)')
ylabel(sp(2),'\itH_s (m)')
xlabel(sp(1),'Model Time (s)')
xlabel(sp(3),'Model Time (s)')
xlabel(sp(2),'Across-shore Distance (m)')
xlabel(sp(4),'Across-shore Distance (m)')

%Save figure
export_fig(ff,[workingDir  'HR_LR_freeSurf_Hs'],'-png','-r600','-nocrop')