% Plot total time-averaged turbulence for different model scenarios input
% into the script with a csv file. Model types are classified by the wave
% period and are plotted as different colors and symbols.
%
% Updates:
% 06/03/21: Load the free surface time series and normalize turbulence by
% incoming wave energy dissipation (WG1 - WG2).
%
% This is version 2 of this script.
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

modelInfo = 'Models_Hs_0_4m_Tp_4-12.csv';
modelDepthScenarios = {'0.5 m';'1 m';'2 m';'3 m';'4 m'};

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Get modelNames from the model directory
HRmodelRaw = dir([HRmodelDir '*rawData_V1*']);
HRmodelRaw = HRmodelRaw(~ismember({HRmodelRaw.name},{'.','..'}));
HRmodelRaw = {HRmodelRaw.name};

HRmodelFree = dir([HRmodelDir '*freeSurf_V1*']);
HRmodelFree = HRmodelFree(~ismember({HRmodelFree.name},{'.','..'}));
HRmodelFree = {HRmodelFree.name};

LRmodelRaw = dir([LRmodelDir '*rawData_V1*']);
LRmodelRaw = LRmodelRaw(~ismember({LRmodelRaw.name},{'.','..'}));
LRmodelRaw = {LRmodelRaw.name};

LRmodelFree = dir([LRmodelDir '*freeSurf_V1*']);
LRmodelFree = LRmodelFree(~ismember({LRmodelFree.name},{'.','..'}));
LRmodelFree = {LRmodelFree.name};

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   600   500]);
pp1 = length(modelInfo.scenarioNumber);
cc = [0.6471 0 0.1490;...
    0.6471 0 0.1490;...
    0.6471 0 0.1490;...
    0.6471 0 0.1490;...
    0.6471 0 0.1490;...
    0.1922 0.2118 0.5843;...
    0.1922 0.2118 0.5843;...
    0.1922 0.2118 0.5843;...
    0.1922 0.2118 0.5843;...
    0.1922 0.2118 0.5843];
markerStyle = {'o';'o';'o';'o';'o';'o';'o';'o';'o';'o'};
disp('Loading HR models')
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(HRmodelRaw','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    if isempty(modelID)
        continue
    else
        fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),HRmodelRaw{modelID})
        load([HRmodelDir HRmodelRaw{modelID}]);
        HR = data;clear data
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
        F1 = trapz(1025*9.81*WG1(2:end).*Cg(2:end));
        F2 = trapz(1025*9.81*WG2(2:end).*Cg(2:end));
        D = (F1-F2)/(63-31); %wave energy dissipation between WG1 and WG2 (Huang et al. 2012)
        
        %Average turbulence in time
        epsAvg = mean(HR.epsilon.epsilon,2);
        x = HR.epsilon.x(:,1);
        z = HR.epsilon.z(:,1);
        
        %     bins = linspace(min(epsAvg),max(epsAvg),1000);
        %     [N,edges,histIDs] = histcounts(epsAvg,bins);
        %     cc = brewermap(length(N),'*RdYlBu');
        %     for j = 1:length(N)
        %         eachBin = find(histIDs == j);
        %         plot3(x(eachBin),z(eachBin),epsAvg(eachBin),'.','color',cc(j,:))
        %         hold on
        %     end
        %     view(0,90)
        
        zIDs = find(z<-0.2); %remove near-surface turbulence in final results
        epsAvg = epsAvg(zIDs);
        
        epsTotal = trapz(epsAvg)/D;
        pp1(i) = plot(modelInfo.h(i),epsTotal,markerStyle{i},'color','k','markerfacecolor',cc(i,:));
        hold on
        
        clear HR
    end
end

pp2 = length(modelInfo.scenarioNumber);
markerStyle = {'d';'d';'d';'d';'d';'d';'d';'d';'d';'d'};
disp('Loading LR models')
%Loop through model scenarios and plot the data
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(LRmodelRaw','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    if isempty(modelID)
        continue
    else
        fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),LRmodelRaw{modelID})
        load([LRmodelDir LRmodelRaw{modelID}]);
        LR = data;clear data
        load([LRmodelDir LRmodelFree{modelID}]);
        FS = data;clear data
        
        %Calculate wave energy dissipation @ WG1 and WG2
        [WG1,~] = pwelch(FS.WG1.eta,[],[],[],8);
        [WG2,F] = pwelch(FS.WG2.eta,[],[],[],8);
        
        k=wavek(F,modelInfo.h(i));
        kh = k*modelInfo.h(i);
        c = sqrt(9.81.*tanh(kh))./sqrt(k);
        n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
        Cg = trapz(n(2:end))*trapz(c(2:end));
        F1 = trapz(WG1(2:end))*Cg*(9.81)*1025;
        F2 = trapz(WG2(2:end))*Cg*(9.81)*1025;
        D = (F1-F2)/(63-31); %wave energy dissipation between WG1 and WG2 (Huang et al. 2012)   
        
        %Average turbulence in time
        epsAvg = mean(LR.epsilon.epsilon,2);
        x = LR.epsilon.x(:,1);
        z = LR.epsilon.z(:,1);
        
        %     bins = linspace(min(epsAvg),max(epsAvg),1000);
        %     [N,edges,histIDs] = histcounts(epsAvg,bins);
        %     cc = brewermap(length(N),'*RdYlBu');
        %     for j = 1:length(N)
        %         eachBin = find(histIDs == j);
        %         plot3(x(eachBin),z(eachBin),epsAvg(eachBin),'.','color',cc(j,:))
        %         hold on
        %     end
        %     view(0,90)
        
        zIDs = find(z<-0.2); %remove near-surface turbulence in final results
        epsAvg = epsAvg(zIDs);
        
        epsTotal = trapz(epsAvg)/D;
        pp2(i) = plot(modelInfo.h(i),epsTotal,markerStyle{i},'color','k','markerfacecolor',cc(i,:));
        hold on
        
        clear LR
    end
end
% Legend
leg = legend([pp1(1) pp1(end) pp2(1) pp2(end)],{'HR Tp = 4 s';'HR Tp = 12 s';'LR Tp = 4 s';'LR Tp = 12 s'});

%Global adjustments
set(gca,'yscale','log',...
    'ylim',[10^-6 10^2],...
    'xtick',0.5)

%Labeling
ylabel('$\overline{\epsilon}/\overline{\epsilon_0} \quad (-)$','interpreter','latex')
xlabel('$\textit{h} \quad (m)$','interpreter','latex')

prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.01 0.01])

%Save figure
% export_fig(ff,[workingDir  'HR_LR_Tp_4s_12s_TotalAvgd_Turbulence'],'-png','-r600','-nocrop')