% Plot bathymetry and water level settings for the HR and LR domains. This
% figure was developed for ICRS 2021 but may also be included in Paper1: 2D
% Wave Flume. This figure shows the 2m depth model as an example. It should
% be noted in the presentation/paper that the other depth models were
% created by raising or lowering the bathymetry while keeping the WSE
% constant at zero.
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
hrBathyPath = 'e:\BKN-FIELD\Models\Paper1_2DwaveFlume\ModelRuns\Calibration\MKK_CombinedHR_V27\MKKextrude\';
lrBathyPath = 'e:\BKN-FIELD\Models\Paper1_2DwaveFlume\ModelRuns\Calibration\MKK_CombinedLR_V2\MKKextrude\';
bathyFileName = 'bathymetry_sample.stl';
saveFigDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\Figures\';

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80 850 500]);
%Define instrument and sample locations
% waveGaugeNames = {'WG1';'WG2';'WG3';'WG4';'WG5';'WG6';'WG7'};
waveGaugeNames = {'WG1';'WG2';'WG3';'WG4';'WG5';'WG6';'WG7';'WG8';'WG9';'WG10';'WG11';'WG12';'WG13';};
HRaqd = 72;
% HRwaveGauges = [31 66 68 70 72 74 76];
HRwaveGauges = [31 65 66 67 68 69 70 71 72 73 74 75 76];
LRaqd = 69;
% LRwaveGauges = [31 63 65 67 69 71 73];
LRwaveGauges = [31 63 64 65 66 67 68 69 70 71 72 73 74];

%HR Domain goes on top
sp(1) = subplot(211);
p(1) = plot3(linspace(29,86,10),ones(1,10)*21.9,zeros(1,10),'color','k','linewidth',1);
hold on
%Load and plot STL
data = stlread([hrBathyPath bathyFileName]);
trimesh(data,'FaceColor','k','EdgeColor','k','linewidth',1.5)
%Plot Aquadopp location
p(2) = plot3(ones(1,10)*HRaqd,repmat(21.9,1,10),linspace(-1.8,0.5,10),'--b','linewidth',1);
p(3) = plot3(HRaqd,21.9,-0.95,'^k','markerfacecolor','k','markersize',6);
%Plot Wave Gauge locations
p(4) = plot3(HRwaveGauges,repmat(21.9,1,length(HRwaveGauges)),repmat(0.5,1,length(HRwaveGauges)),'+k','markersize',8);
p(5) = plot3(HRwaveGauges,repmat(21.9,1,length(HRwaveGauges)),repmat(0.5,1,length(HRwaveGauges)),'ok','markersize',4);
view(0,0)
set(gca,'zlim',[-3 2])

%LR Domain goes on bottom
sp(2) = subplot(212);
p(6) = plot3(linspace(29,86,10),ones(1,10)*21.9,zeros(1,10),'color','k','linewidth',1);
hold on
%Load and plot STL
data = stlread([lrBathyPath bathyFileName]);
trimesh(data,'FaceColor','k','EdgeColor','k','linewidth',1.5)
%Plot Aquadopp location
p(7) = plot3(ones(1,10)*LRaqd,repmat(10.9,1,10),linspace(-1.8,0.5,10),'--b','linewidth',1);
p(8) = plot3(LRaqd,10.9,-0.95,'^k','markerfacecolor','k','markersize',6);
%Plot Wave Gauge locations
p(9) = plot3(LRwaveGauges,repmat(10.9,1,length(LRwaveGauges)),repmat(0.5,1,length(LRwaveGauges)),'+k','markersize',8);
p(10) = plot3(LRwaveGauges,repmat(10.9,1,length(LRwaveGauges)),repmat(0.5,1,length(LRwaveGauges)),'ok','markersize',4);
view(0,0)
set(gca,'zlim',[-3 2])

%Figure annotations
plot3(sp(1),[65 76 76 65 65], [21.9 21.9 21.9 21.9 21.9], [-3 -3 2 2 -3],'--r','linewidth',1)
plot3(sp(2),[63 74 74 63 63], [21.9 21.9 21.9 21.9 21.9], [-3 -3 2 2 -3],'--r','linewidth',1)

%Global plot adjustments
set(sp,'xlim',[30 85])
set(sp(1),'xticklabel',[])
set(sp(2),'xticklabel',{'0';'5';'10';'15';'20';'25';'30';'35';'40';'45';'50';'55'})

%Labeling
zlabel(sp(1),'z (m)')
zlabel(sp(2),'z (m)')
xlabel(sp(2),'Across-shore Distance (m)')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.005 0.005])

%Text annotations
text(sp(1),31,21.9,1.5,'High Relief','fontsize',11)
text(sp(1),HRaqd+1,21.9,-0.95,'HR Aquadopp','fontsize',8)
% text(sp(1),HRwaveGauges-0.8,repmat(21.9,1,length(HRwaveGauges)),repmat(0.9,1,length(HRwaveGauges)),waveGaugeNames,'fontsize',6)
text(sp(2),31,10.9,1.5,'Low Relief','fontsize',11)
text(sp(2),LRaqd+1,10.9,-0.95,'LR Aquadopp','fontsize',8)
% text(sp(2),LRwaveGauges-0.8,repmat(10.9,1,length(LRwaveGauges)),repmat(0.9,1,length(LRwaveGauges)),waveGaugeNames,'fontsize',6)


set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(ff,[saveFigDir  'HR_LR_calibration_domains_with_annotations_V4'],'-dpdf','-r0')
%export_fig(ff,[saveFigDir  'HR_LR_calibration_domains_with_annotations_V3'],'-pdf','-nocrop','-nofontswap')
