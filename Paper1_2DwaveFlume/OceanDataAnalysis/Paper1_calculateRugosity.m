%Calculate rugosity of the HR and LR patches
%
% Updates:
% 10/8/21: added sigma_r and k_w calcs
% This is version 1 of this script.
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

%% HR
data = stlread([hrBathyPath bathyFileName]);
[bathyX,sortID] = sort(data.Points(:,1));
bathyY = data.Points(sortID,3);
% HRpatch = [64 76];
HRpatch = [65 71];
idx = find(bathyX>=HRpatch(1)&bathyX<=HRpatch(2));
xy = [bathyX(idx) bathyY(idx)];
Lchain = sum(sqrt(sum(abs(diff(xy)).^2,2)));
Dchain = sum(sqrt(sum(abs(diff(HRpatch)).^2,2))); %this should equal HRpatch(2)-HRpatch(1)
Rugosity = Lchain/Dchain;
fprintf('HR rugosity index: %0.2f\n',Rugosity)
sigmar = std(bathyY(idx));
kw = 4*sigmar;
fprintf('sigma_r = %0.3f, k_w = %0.3f\n',sigmar,kw)
%% LR
data = stlread([lrBathyPath bathyFileName]);
[bathyX,sortID] = sort(data.Points(:,1));
bathyY = data.Points(sortID,3);
% LRpatch = [62 74];
LRpatch = [69.5 74];
idx = find(bathyX>=LRpatch(1)&bathyX<=LRpatch(2));
xy = [bathyX(idx) bathyY(idx)];
Lchain = sum(sqrt(sum(abs(diff(xy)).^2,2)));
Dchain = sum(sqrt(sum(abs(diff(LRpatch)).^2,2))); %this should equal LRpatch(2)-LRpatch(1)
Rugosity = Lchain/Dchain;
fprintf('LR rugosity index: %0.2f\n',Rugosity)
sigmar = std(bathyY(idx));
kw = 4*sigmar;
fprintf('sigma_r = %0.3f, k_w = %0.3f\n',sigmar,kw)
