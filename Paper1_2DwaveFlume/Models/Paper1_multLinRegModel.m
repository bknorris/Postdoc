% Run a multiple linear regression model on the HR and LR datasets to
% determine which forcing condition has the greatest effect on turbulence.
%
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
modelInfo = 'Models_allHR_LR.csv';

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%% Do the HR models
HR = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V4.mat']);

%Extract the epsilon data
patchLo = 35;
patchHi = 41;
zPatch = [-0.9 -1.9 -2.9 -3.9];
upper = 25;
eps = zeros(length(HR.Hs),1);
for j = 1:length(HR.scenarioNumber)
    eps1 = squeeze(HR.eps(:,:,j));
    xBins = squeeze(HR.xBins(:,:,j));
    zBins = squeeze(HR.zBins(:,:,j));
    [~,xLo] = min(abs(xBins-patchLo));
    [~,xHi] = min(abs(xBins-patchHi));
    if HR.h(j) == 1
        [~,zHi] = min(abs(zBins-zPatch(1)));
    elseif HR.h(j) == 2
        [~,zHi] = min(abs(zBins-zPatch(2)));
    elseif HR.h(j) == 3
        [~,zHi] = min(abs(zBins-zPatch(3)));
    elseif HR.h(j) == 4
        [~,zHi] = min(abs(zBins-zPatch(4)));
    end
    
    %Find the "bed" as the last value in the column that is a nan
    epscrop = eps1(xLo:xHi,:);[m,n] = size(epscrop);
    epsfilt = NaN(m,n);
    for k = 1:m
        nanID = find(isnan(epscrop(k,1:end-1)),1,'last');
        epsfilt(k,nanID+1:nanID+upper) = epscrop(k,nanID+1:nanID+upper);
    end
    eps(j) = nanmean(nanmean(epsfilt));
end

%Run stepwiseLM
disp('Running MLR on HR data')
vals = {'h';'Hs';'Tp';'eps'};
Hs = HR.Hs';
Tp = HR.Tp';
h = HR.h';
data = [ones(size(eps)) h Hs Tp];
y = eps;
[b,bint,r,rint,stats]  = regress(zscore(y),[data(:,1) zscore(data(:,2:end))],0.05);
fprintf('R-squared: %0.3f\n',stats(1))
lma = stepwiselm(zscore(data(:,2:end)),zscore(y),'constant','Upper','linear',...
    'penter',0.05,'premove',0.1,'varnames',vals)
fprintf('\n\n')
%% Do the LR models
LR = load([LRmodelDir 'LR_postProcess_epsilon_allHR_LR_V4.mat']);

%Extract the epsilon data
patchLo = 35;
patchHi = 41;
zPatch = [-0.9 -1.9 -2.9 -3.9];
upper = 25;
eps = zeros(length(LR.Hs),1);
for j = 1:length(LR.scenarioNumber)
    eps1 = squeeze(LR.eps(:,:,j));
    xBins = squeeze(LR.xBins(:,:,j));
    zBins = squeeze(LR.zBins(:,:,j));
    [~,xLo] = min(abs(xBins-patchLo));
    [~,xHi] = min(abs(xBins-patchHi));
    if LR.h(j) == 1
        [~,zHi] = min(abs(zBins-zPatch(1)));
    elseif LR.h(j) == 2
        [~,zHi] = min(abs(zBins-zPatch(2)));
    elseif LR.h(j) == 3
        [~,zHi] = min(abs(zBins-zPatch(3)));
    elseif LR.h(j) == 4
        [~,zHi] = min(abs(zBins-zPatch(4)));
    end
    
    %Find the "bed" as the last value in the column that is a nan
    epscrop = eps1(xLo:xHi,:);[m,n] = size(epscrop);
    epsfilt = NaN(m,n);
    for k = 1:m
        nanID = find(isnan(epscrop(k,1:end-1)),1,'last');
        epsfilt(k,nanID+1:nanID+upper) = epscrop(k,nanID+1:nanID+upper);
    end
    eps(j) = nanmean(nanmean(epsfilt));
end

%Run stepwiseLM
disp('Running MLR on LR data')
vals = {'h';'Hs';'Tp';'eps'};
Hs = LR.Hs';
Tp = LR.Tp';
h = LR.h';
data = [ones(size(eps)) h Hs Tp];
y = eps;
[b,bint,r,rint,stats]  = regress(zscore(y),[data(:,1) zscore(data(:,2:end))],0.01);
fprintf('R-squared: %0.3f\n',stats(1))
lma = stepwiselm(zscore(data(:,2:end)),zscore(y),'constant','Upper','linear',...
    'penter',0.05,'premove',0.1,'varnames',vals)

