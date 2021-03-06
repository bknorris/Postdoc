% Create a plot of normalized TKE dissipation for HR vs LR models to
% directly compare the two.
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

%Load processed model data
HR = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V4.mat']);
LR = load([LRmodelDir 'LR_postProcess_epsilon_allHR_LR_V4.mat']);

%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   500   400]);
wavePeriods = unique(modelInfo.Tp);
lineStyle = {'-';'--';'-.';':';'-'};
markers = {'o';'s';'d';'^';'p'};
sp = zeros(2,1);
cc = hex2rgb({'#1e5cb3';'#04a1e6';'#ffdf28';'#f58b35';'#bf171f'}); %blue-yellow-red colormap based on GMT panopoly
xLin = linspace(10E-4,10E0,10);
yLin = linspace(10E-4,10E0,10);
plot(xLin,yLin,'--k','linewidth',1.5);hold on

for i = 1:length(wavePeriods)
    %Do the HR models first; make sure these are the same as
    %Paper1_plotDepthProfiles.m!
    patchLo = 35;
    patchHi = 41;
    zPatch = [-0.9 -1.9 -2.9 -3.9];
    
    %Organize plots by water depth!
    modelID = find(HR.Tp == wavePeriods(i));
    HRepsNorm = zeros(length(modelID),1);
    
    for j = 1:length(modelID)
        eps = squeeze(HR.eps(:,:,modelID(j)));
        Uw = squeeze(HR.Umag(:,:,modelID(j)));
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
          HRepsNorm(j) = nanmean(nanmean(eps(xLo:xHi,:)./(Uw(xLo:xHi,:).^3)));
    end
    
    %Then do LR models
    patchLo = 40.5;
    patchHi = 44;
    zPatch = [-0.92 -1.92 -2.92 -3.92];
    
    %Organize plots by water depth!
    modelID = find(LR.Tp == wavePeriods(i));
    LRepsNorm = zeros(length(modelID),1);
    for j = 1:length(modelID)
        eps = squeeze(LR.eps(:,:,modelID(j)));
        Uw = squeeze(LR.Umag(:,:,modelID(j)));
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
          LRepsNorm(j) = nanmean(nanmean(eps(xLo:xHi,:)./(Uw(xLo:xHi,:).^3)));
    end

    x = HRepsNorm;
    y = LRepsNorm;
    plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
%     B = polyfit(x,log10(y),1);
%     Yfit = 10.^(B(1)*linspace(0.1,0.7,length(x))+B(2));
%     plot(linspace(0.1,0.7,length(x)),Yfit,'-.',...
%         'color',cc(i,:),'linewidth',1)

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
    'xlim',[5E-3 10E-1],'ylim',[5E-3 10E-1])
set(leg,'position',[0.18 0.7 0.1 0.2])

%Labeling
xlabel('High Relief $\langle \overline{\varepsilon}\rangle/\langle \overline{{U_w}^3} \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
ylabel('Low Relief $\langle \overline{\varepsilon}\rangle/\langle \overline{{U_w}^3} \rangle \  \mathrm{(m^{-1})}$','interpreter','latex')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.005 0.005])

set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% print(ff,[workingDir 'AllModels_NormalizedDissipation_HRvsLR'],'-dpdf','-r0')