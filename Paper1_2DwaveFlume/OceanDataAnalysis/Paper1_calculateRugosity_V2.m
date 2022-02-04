%Calculate rugosity of the HR and LR patches
%
% Updates:
% 10/8/21:  added sigma_r and k_w calcs
% 10/25/21: added height and wavelength calcs and plot routine to plot
%           bathy profile on the same figure
%
% This is version 2 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%Define data paths
hrBathyPath = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\Calibration\MKK_CombinedHR_V27\MKKextrude\';
lrBathyPath = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\Calibration\MKK_CombinedLR_V2\MKKextrude\';
bathyFileName = 'bathymetry_sample.stl';
saveFigDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\Figures\';
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80 500 250]);
cc = [237 28 36;29 140 204]./255;
%% HR
data = stlread([hrBathyPath bathyFileName]);
[bathyX,sortID] = sort(data.Points(:,1));
bathyY = data.Points(sortID,3);
HRpatch = [67 71.2];
idx = find(bathyX>=HRpatch(1)&bathyX<=HRpatch(2));
xy = [bathyX(idx) bathyY(idx)];
Lchain = sum(sqrt(sum(abs(diff(xy)).^2,2)));
Dchain = sum(sqrt(sum(abs(diff(HRpatch)).^2,2))); %this should equal HRpatch(2)-HRpatch(1)
Rugosity = Lchain/Dchain;
fprintf('HR rugosity index: %0.2f\n',Rugosity)
bathyY = bathyY(idx);bathyX = bathyX(idx); %redefine bathy
sigmar = std(bathyY);
kw = 4*sigmar;
fprintf('sigma_r = %0.3f, k_w = %0.3f\n',sigmar,kw)
%Get mean wavelength and max height of roughness
% %Use findpeaks
[~,maxLoc] = findpeaks(bathyY,'MinPeakProminence',0.005);
[~,minLoc] = findpeaks(-1*bathyY,'MinPeakProminence',0.005);
minLoc = [1; minLoc]; %include first value because it is a local minimum
minima = bathyY(minLoc);maxima = bathyY(maxLoc);
heightDiffs = maxima-minima;
sigHeight = mean(maxk(heightDiffs,round(length(heightDiffs)/3,0)));
minima = bathyX(minLoc);maxima = bathyX(maxLoc);
lengthDiffs = maxima-minima;
sigLength = mean(maxk(lengthDiffs,round(length(heightDiffs)/3,0)));
fprintf('Roughness length = %0.3f m, Roughness height = %0.3f m\n',sigLength,sigHeight)    

% figure
% plot(bathyX,bathyY,'k')
% hold on
% plot(bathyX(maxLoc),bathyY(maxLoc),'^r')
% plot(bathyX(minLoc),bathyY(minLoc),'ob')
% legend
% xlabel('Across-shore (m)')
% ylabel('Elevation (m)')

%Plot routine
xs = linspace(0,4.2,length(bathyY));
p(1) = plot(xs,smooth(bathyY,3),'color',cc(1,:),'linewidth',1.2); hold on

%% LR
data = stlread([lrBathyPath bathyFileName]);
[bathyX,sortID] = sort(data.Points(:,1));
bathyY = data.Points(sortID,3);
LRpatch = [69 73.2];
idx = find(bathyX>=LRpatch(1)&bathyX<=LRpatch(2));
xy = [bathyX(idx) bathyY(idx)];
Lchain = sum(sqrt(sum(abs(diff(xy)).^2,2)));
Dchain = sum(sqrt(sum(abs(diff(LRpatch)).^2,2))); %this should equal LRpatch(2)-LRpatch(1)
Rugosity = Lchain/Dchain;
fprintf('LR rugosity index: %0.2f\n',Rugosity)
bathyY = bathyY(idx);bathyX = bathyX(idx); %redefine bathy
sigmar = std(bathyY);
kw = 4*sigmar;
fprintf('sigma_r = %0.3f, k_w = %0.3f\n',sigmar,kw)
%Get mean wavelength and max height of roughness
%Use findpeaks
[~,maxLoc] = findpeaks(bathyY,'MinPeakProminence',0.005);
[~,minLoc] = findpeaks(-1*bathyY,'MinPeakProminence',0.005);
minLoc = [1; minLoc]; %include first value because it is a local minimum
maxLoc = [maxLoc; length(bathyY)]; %include last value because it is a local maximum
minima = bathyY(minLoc);maxima = bathyY(maxLoc);
heightDiffs = maxima-minima;
sigHeight = mean(maxk(heightDiffs,round(length(heightDiffs)/3,0)));
minima = bathyX(minLoc);maxima = bathyX(maxLoc);
lengthDiffs = maxima-minima;
sigLength = mean(maxk(lengthDiffs,round(length(heightDiffs)/3,0)));
fprintf('Roughness length = %0.3f m, Roughness height = %0.3f m\n',sigLength,sigHeight)    

figure
plot(bathyX,bathyY,'k')
hold on
plot(bathyX(maxLoc),bathyY(maxLoc),'^r')
plot(bathyX(minLoc),bathyY(minLoc),'ob')
legend
xlabel('Across-shore (m)')
ylabel('Elevation (m)')

%Plot routine
xs = linspace(0,4.2,length(bathyY));
p(2) = plot(xs,smooth(bathyY,3),'color',cc(2,:),'linewidth',1.2);
leg = legend(p,{'High Relief';'Low Relief'},'location','northeast');
set(leg,'box','off')

set(gca,'xlim',[0 4.2])
xlabel('Across-shore Distance (m)')
ylabel('z (m)')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.008 0.008])
set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(ff,[saveFigDir 'PatchBathyProfiles'],'-dpdf','-r0')

