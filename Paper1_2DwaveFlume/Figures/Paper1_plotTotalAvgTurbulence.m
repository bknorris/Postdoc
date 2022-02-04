% Plot total time-averaged turbulence for different model scenarios input
% into the script with a csv file. Model types are classified by the wave
% period and are plotted as different colors and symbols.
%
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
HRmodelRaw = dir([HRmodelDir '*rawData*']);
HRmodelRaw = HRmodelRaw(~ismember({HRmodelRaw.name},{'.','..'}));
HRmodelRaw = {HRmodelRaw.name};

LRmodelRaw = dir([LRmodelDir '*rawData_V1*']);
LRmodelRaw = LRmodelRaw(~ismember({LRmodelRaw.name},{'.','..'}));
LRmodelRaw = {LRmodelRaw.name};

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
        
        epsTotal = trapz(epsAvg);
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
        
        epsTotal = trapz(epsAvg);
        pp2(i) = plot(modelInfo.h(i),epsTotal,markerStyle{i},'color','k','markerfacecolor',cc(i,:));
        hold on
        
        clear LR
    end
end
%
% %Create legend
% cb = colorbar('location','southoutside');
% set(cb,'position',[0.16 0.08 0.7 0.02],...
%     'tickdir','out',...
%     'ticklength',0.01,...
%     'linewidth',1.5)
%
% % %Global Adjustments
% set([sp1(1) sp2(1)],'ylim',[-1.1 -0.2],'ytick',-1:0.4:-0.2)
% set([sp1(2) sp2(2)],'ylim',[-1.6 -0.4])
% set([sp1(3) sp2(3)],'ylim',[-2.6 -1])
% set([sp1(4) sp2(4)],'ylim',[-3.6 -2])
% set([sp1(5) sp2(5)],'ylim',[-4.6 -3])
%
% set(sp1(1),'position',[0.1 0.8 0.38 0.12],'xticklabel',[])
% set(sp1(2),'position',[0.1 0.65 0.38 0.12],'xticklabel',[])
% set(sp1(3),'position',[0.1 0.5 0.38 0.12],'xticklabel',[])
% set(sp1(4),'position',[0.1 0.35 0.38 0.12],'xticklabel',[])
% set(sp1(5),'position',[0.1 0.2 0.38 0.12])
%
% set(sp2(1),'position',[0.55 0.8 0.38 0.12],'xticklabel',[],'yticklabel',[])
% set(sp2(2),'position',[0.55 0.65 0.38 0.12],'xticklabel',[],'yticklabel',[])
% set(sp2(3),'position',[0.55 0.5 0.38 0.12],'xticklabel',[],'yticklabel',[])
% set(sp2(4),'position',[0.55 0.35 0.38 0.12],'xticklabel',[],'yticklabel',[])
% set(sp2(5),'position',[0.55 0.2 0.38 0.12],'yticklabel',[])
%
% % %Labeling
% title(sp1(1),'HR Models')
% title(sp2(1),'LR Models')
% ylabel(sp1(1),'z (m)')
% ylabel(sp1(2),'z (m)')
% ylabel(sp1(3),'z (m)')
% ylabel(sp1(4),'z (m)')
% ylabel(sp1(5),'z (m)')
% xlabel(sp1(5),'Across-shore Distance (m)')
% xlabel(sp2(5),'Across-shore Distance (m)')
% xlabel(cb,'$\overline{\epsilon} \quad (m^2/s^3)$','interpreter','latex')
% %
% prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.01 0.01])
%
% %Save figure
% export_fig(ff,[workingDir  'HR_LR_Tp_4s_TimeAvgd_Turbulence'],'-png','-r600','-nocrop')