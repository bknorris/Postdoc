% Plot depth-averaged velocities from the calibration "Mesh size
% refinement" models. 
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all
dataPath = 'c:\Users\user\Documents\Models\Paper2_OptimizingRestoration\ModelRuns\MeshScaleTests\DataAnalysis\';
figPath = 'c:\Users\user\Documents\Models\Paper2_OptimizingRestoration\Figures\'; %where to save figs
modelLegend = {'Level 1 [56 - 28 mm]';'Level 2 [56 - 14 mm]';'Level 3 [56 - 7 mm]'};
modelFolders = dir([dataPath '*rawData*']);
modelFolders = {modelFolders.name};
Zbins = linspace(-0.028,0,40);


%Plot routine
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   500   400]);
hold on;
p = zeros(1,length(modelFolders));
cc = brewermap(length(modelFolders),'Blues');
for i = 1:length(modelFolders)
    data = load([dataPath modelFolders{i}]);  
    [~,dataLength] = size(data.U.x);
    zBin = zeros(32,dataLength);
    UmagBin = zeros(32,dataLength);
    %Load the data
    for j = 1:dataLength
        z = data.U.z(:,j);
        Umag = sqrt(data.U.Ux(:,j).^2);
        %bin the epsilon data by the z axis
        [~,~,loc] = histcounts(z,Zbins);
        lz = find(loc>0);
        zBin(:,j) = accumarray(loc(lz),z(lz))./accumarray(loc(lz),1);
        UmagBin(:,j) = accumarray(loc(lz),Umag(lz))./accumarray(loc(lz),1);
    end
    p(i) = plot(smooth(nanmean(UmagBin,2),8),nanmean(zBin,2).*-1000,...
        'linewidth',1.5,'color',cc(i,:));   
end
set(gca,'xlim',[0.01 0.04],'ydir','reverse','ylim',[0 28])

%Figure annotations
leg = legend(p,modelLegend,'location','northwest');
 
% %Global plot adjustments
% set(sp,'xlim',[30 85])
% set(sp(1),'xticklabel',[])
% set(sp(2),'xticklabel',{'0';'5';'10';'15';'20';'25';'30';'35';'40';'45';'50';'55'})

% %Labeling
ylabel('Depth (mm)','interpreter','latex')
xlabel('$|\overline{U}|$ (m/s)','interpreter','latex')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.005 0.005])
 
set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(ff,[figPath  'MeshSize_calibration'],'-dpng','-r0')
