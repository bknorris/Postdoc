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

%Load processed model data
HR = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V6.mat']);
LR = load([LRmodelDir 'LR_postProcess_epsilon_allHR_LR_V6.mat']);

ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   500   400]);
cc = flipud(hex2rgb({'#032760';'#0057d2';'#57b0ff';'#73dcff'})); %A dark -> light blue colormap (based on Bath 112)
% cc = [32 80 255;134 217 255;255 196 0;213 0 0]./255; %A red orange light blue blue colormap (based on GMT Panopoly)
markers = {'^';'s';'o';'d'};
lineStyle = {'-';'--';'-.';':'};
% sp(1) = subplot(121);
xs = repmat(linspace(-0.6,0.6,4)',1,5);
waterDepths = unique(modelInfo.h);
wavePeriods = unique(modelInfo.Tp);
upper = 15;
for i = 1:length(waterDepths)
    %Do the HR models first0
    patchLo = 35;
    patchHi = 41;
    zPatch = [-0.9 -1.9 -2.9 -3.9];
    modelID = find(HR.h == waterDepths(i));

    %EXPERIMENTAL: try filtering results by their proximity to the bed
    for j = 1:length(modelID)
        Uw1 = squeeze(HR.Uw(:,:,modelID(j)));
        xBins = squeeze(HR.xBins(:,:,modelID(j)));
        zBins = squeeze(HR.zBins(:,:,modelID(j)));
        [~,xLo] = min(abs(xBins-patchLo));
        [~,xHi] = min(abs(xBins-patchHi));
        if HR.h(modelID(j)) == 1
            [~,zHi] = min(abs(zBins-zPatch(1)));
        elseif HR.h(modelID(j)) == 2
            [~,zHi] = min(abs(zBins-zPatch(2)));
        elseif HR.h(modelID(j)) == 3
            [~,zHi] = min(abs(zBins-zPatch(3)));
        elseif HR.h(modelID(j)) == 4
            [~,zHi] = min(abs(zBins-zPatch(4)));
        end
        
        %Find the "bed" as the last value in the column that is a nan
        uwcrop = Uw1(xLo:xHi,:);[m,n] = size(uwcrop);
        uwfilt = NaN(m,n);
        for k = 1:m
            nanID = find(isnan(uwcrop(k,1:end-1)),1,'last');
            uwfilt(k,nanID+1:nanID+upper) = uwcrop(k,nanID+1:nanID+upper);
        end
        uw(j) = nanmean(nanmean(uwfilt));
        Dmean(j) = mean(HR.D(:,modelID(j)),1);
        Tp(j) = HR.Tp(modelID(j));
    end
    fe = abs(Dmean)./((1025*(2/3*pi)*uw.^3));
    
    x = Tp;
    y = fe;
%     y = HsNorm;
%     plot(x,y,sprintf('%s',markers{i}),...
%         'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
        edges = linspace(0,24,6);
        bins = discretize(x,edges);
        [meanVal,maxVal,stDev] = grpstats(y,bins,{@mean, @max, @std});
        errorbar(wavePeriods'+xs(i,:),meanVal,stDev,...
            sprintf('%s%s',markers{i},lineStyle{i}),'color',cc(i,:),...
            'markerfacecolor',cc(i,:),'markersize',6,'linewidth',1.2,'capsize',0);
    hold on;
end
