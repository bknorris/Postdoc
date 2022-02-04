% Plot depth-averaged velocities from the calibration "Mesh size
% refinement" models. 
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all
dataPath = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\Calibration\Mesh_scale_vs_turbulence\';
figPath = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\Figures\'; %where to save figs
modelLegend = {'Level 1 [0.3 - 0.15 m]';'Level 2 [0.3 - 0.075 m]';'Level 3 [0.3 - 0.0375 m]';...
    'Level 4 [0.3 - 0.0187 m]';'Level 5 [0.3 - 0.01 m]'};
modelFolders = dir([dataPath '*VelocityData*']);
modelFolders = {modelFolders.name};

%Plot routine
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   500   400]);
hold on;
p = zeros(1,length(modelFolders));
cc = brewermap(length(modelFolders),'Blues');
for i = 1:length(modelFolders)
     load([dataPath modelFolders{i}])   
     p(i) = plot(smooth(nanmean(data.UmagBin,2),8),nanmean(data.zBin,2).*-1,...
        'linewidth',1.5,'color',cc(i,:));
end
set(gca,'ydir','reverse','ylim',[0.8 2])

%Figure annotations
leg = legend(p,modelLegend,'location','northwest');
 
% %Global plot adjustments
% set(sp,'xlim',[30 85])
% set(sp(1),'xticklabel',[])
% set(sp(2),'xticklabel',{'0';'5';'10';'15';'20';'25';'30';'35';'40';'45';'50';'55'})

% %Labeling
ylabel('Depth (m)','interpreter','latex')
xlabel('$|\overline{U}|$ (m/s)','interpreter','latex')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.005 0.005])
 
set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(ff,[figPath  'MeshSize_calibration'],'-dpdf','-r0')
