% Compare the populations of HR and LR models with the same wave period to
% see if they are statistically different (t-test).
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
wavePeriods = unique(modelInfo.Tp);

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   500   400]);
lineStyle = {'-';'--';'-.';':';'-'};
markers = {'o';'s';'d';'^';'p'};
sp = zeros(2,1);
cc = hex2rgb({'#1e5cb3';'#04a1e6';'#ffdf28';'#f58b35';'#bf171f'}); %blue-yellow-red colormap based on GMT panopoly
xLin = linspace(10E-4,10E0,10);
yLin = linspace(10E-4,10E0,10);
plot(xLin,yLin,'--k','linewidth',1.5);hold on

%Load processed model data
HR = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V4.mat']);
LR = load([LRmodelDir 'LR_postProcess_epsilon_allHR_LR_V4.mat']);
upper = 45;
for i = 1:length(wavePeriods)
    
    %Do the HR models first0
    patchLo = 35;
    patchHi = 41;
    zPatch = [-0.9 -1.9 -2.9 -3.9];
    modelID = find(HR.Tp == wavePeriods(i));
    HRepsNorm = zeros(length(modelID),1);
    
    %EXPERIMENTAL: try filtering results by their proximity to the bed
    for j = 1:length(modelID)
        eps1 = squeeze(HR.eps(:,:,modelID(j)));
        Uw1 = squeeze(HR.Umag(:,:,modelID(j)));
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
        epscrop = eps1(xLo:xHi,:);[m,n] = size(epscrop);
        uwcrop = Uw1(xLo:xHi,:);
        epsfilt = NaN(m,n);
        uwfilt = NaN(m,n);
        for k = 1:m
            nanID = find(isnan(epscrop(k,1:end-1)),1,'last');
            epsfilt(k,nanID+1:nanID+upper) = epscrop(k,nanID+1:nanID+upper);
            nanID = find(isnan(uwcrop(k,1:end-1)),1,'last');
            uwfilt(k,nanID+1:nanID+upper) = uwcrop(k,nanID+1:nanID+upper);
        end
        HRepsNorm(j) = nanmean(nanmean(epsfilt./(uwfilt.^3)));
    end
    
    %Then do the LR models
    patchLo = 40.5;
    patchHi = 43.5;
    zPatch = [-0.92 -1.92 -2.92 -3.92];
    modelID = find(LR.Tp == wavePeriods(i));
    LRepsNorm = zeros(length(modelID),1);
    
    %EXPERIMENTAL: try filtering results by their proximity to the bed
    for j = 1:length(modelID)
        eps1 = squeeze(LR.eps(:,:,modelID(j)));
        Uw1 = squeeze(LR.Umag(:,:,modelID(j)));
        xBins = squeeze(LR.xBins(:,:,modelID(j)));
        zBins = squeeze(LR.zBins(:,:,modelID(j)));
        [~,xLo] = min(abs(xBins-patchLo));
        [~,xHi] = min(abs(xBins-patchHi));
        if LR.h(modelID(j)) == 1
            [~,zHi] = min(abs(zBins-zPatch(1)));
        elseif LR.h(modelID(j)) == 2
            [~,zHi] = min(abs(zBins-zPatch(2)));
        elseif LR.h(modelID(j)) == 3
            [~,zHi] = min(abs(zBins-zPatch(3)));
        elseif LR.h(modelID(j)) == 4
            [~,zHi] = min(abs(zBins-zPatch(4)));
        end
        
        %Find the "bed" as the last value in the column that is a nan
        epscrop = eps1(xLo:xHi,:);[m,n] = size(epscrop);
        uwcrop = Uw1(xLo:xHi,:);
        epsfilt = NaN(m,n);
        uwfilt = NaN(m,n);
        for k = 1:m
            nanID = find(isnan(epscrop(k,1:end-1)),1,'last');
            epsfilt(k,nanID+1:nanID+upper) = epscrop(k,nanID+1:nanID+upper);
            nanID = find(isnan(uwcrop(k,1:end-1)),1,'last');
            uwfilt(k,nanID+1:nanID+upper) = uwcrop(k,nanID+1:nanID+upper);
        end
        LRepsNorm(j) = nanmean(nanmean(epsfilt./(uwfilt.^3)));
    end
    
    %Run the two-sided t-test on each sample population
    [h,p,ci,stats] = ttest2(HRepsNorm,LRepsNorm,'vartype','unequal');
    fprintf('Wave Period: %0.0f seconds\n',wavePeriods(i))
    if h == 1
        fprintf('The sample population from HR and LR models are statistically different at the 5%% level\n')
    else
        fprintf('The sample population from HR and LR models are not statistically different at the 5%% level\n')
    end
    fprintf('p-val: %0.2d\n',p)
    fprintf('t-stat: %0.2f\n',stats.tstat)
    fprintf('df: %0.2f\n',stats.df)
    
    %Plot results
    x = HRepsNorm;
    y = LRepsNorm;
    plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
end

%Make some dummy data for the legends
xs = repmat(linspace(1E5,1E6,2),length(wavePeriods),1);
ys = repmat(linspace(min(x),max(x),2),length(wavePeriods),1);
pp3 = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    pp3(i) = plot(xs(i,:),ys(i,:),sprintf('%s',markers{i}),...
        'markerfacecolor',cc(i,:),'markersize',4,'color',cc(i,:),'linewidth',1);hold on;
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end
leg = legend(pp3,legText,'location','southeast','box','off','interpreter','latex');

%Axes and plot adjustments
set(gca,'xscale','log','yscale','log',...
    'xlim',[8E-3 2],'ylim',[8E-3 2])
set(leg,'position',[0.18 0.7 0.1 0.2])

%Labeling
xlabel('High Relief $\langle \overline{\varepsilon}\rangle/\langle \overline{{U_w}^3} \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
ylabel('Low Relief $\langle \overline{\varepsilon}\rangle/\langle \overline{{U_w}^3} \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.005 0.005])

set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(ff,[workingDir 'AllModels_NormalizedDissipation_HRvsLR'],'-dpdf','-r0')

