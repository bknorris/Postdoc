% Create a plot of wave dissipation factor fe versus wave period for
% four different water depths -- This method uses Gon et al., 2020 to
% calculate fe. Currently it works but the results do not look right. 
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Load the data
%Define data paths
workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\Figures\';
HRmodelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\HR_Domain\';
LRmodelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\LR_Domain\';
HRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\Models\';
LRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Models\';

%% Load processed model data -- HR
HRmodelFree = dir([HRmodelDir '*freeSurf_V2*']);
HRmodelFree = HRmodelFree(~ismember({HRmodelFree.name},{'.','..'}));
HRmodelFree = {HRmodelFree.name};

HRmodelBathy = dir([HRbathyDir '*_bathy_0.stl']);
HRmodelBathy = HRmodelBathy(~ismember({HRmodelBathy.name},{'.','..'}));
HRmodelBathy = {HRmodelBathy.name};

modelNames = 'Models_allHR.csv';
fid = fopen([workingDir modelNames]);
header = textscan(fgetl(fid),repmat('%s',1,9),'delimiter',',');
data = textscan(fid,repmat('%f',1,9),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);clear data

disp('Loading HR models')
%Loop through model scenarios and load the data into a structure
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(HRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    FS = load([HRmodelDir HRmodelFree{modelID}]);
    fn = fieldnames(FS);
    h = modelInfo.h(i);
    HR.h(i) = h;
    HR.Tp(i) = modelInfo.Tp(i);
    F = zeros(12,64);
    ubj = zeros(12,64); %orbital velocity from linear wave theory
    ub = zeros(12,1);
    
    %Get local h with bathy file
    whichBathy = find(contains(HRmodelBathy,sprintf('%dm',h)));
    bathy = stlread([HRbathyDir HRmodelBathy{whichBathy}]);
    [bathyX,sortID] = sort(bathy.Points(:,1));
    bathyY = bathy.Points(sortID,3);
    
    %Crop time series based on start/stop times in modelInfo
    times = FS.time;
    timeID = find(times >= modelInfo.tStart(i) & times <= modelInfo.tStop(i));
    
    %Update 09/30/21: Calculation of wave dissipation factor fe from
    %surface elevation spectra (Lowe et al. 2005)
    WGs = FS.waveGauges;
    WGorder = [1 6 7 8 9 10 11 12 13 2 3 4 5]; %wave gauges are out of numerical order in FS
    for j = 1:length(WGorder)
        eta = FS.(fn{WGorder(j)}).eta(timeID);
        eta = repmat(eta,1,5); %make the time series longer
        [SSeta,f] = pwelch(eta,256,[],[],8);
        [~,bathyID] = min(abs(bathyX-WGs(j)));
        localh = -1*bathyY(bathyID); %local h at wave gauge
        
        k = wavek(f,localh);
        kh = k*localh;
        c = sqrt(9.81.*tanh(kh))./sqrt(k);
        n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
        Cg = n.*c;
        cutoff = find(f>=2,1,'first');
        F(j,:) = 1025*9.81*SSeta(2:cutoff).*Cg(2:cutoff);
        
        %Estimate near-bottom orbital velocity with linear wave theory
        aj = sqrt(2*SSeta(2:cutoff));
        omegaj = 2*pi*f(2:cutoff);
        ubj(j,:) = (aj.*omegaj)./sinh(kh(2:cutoff));  
    end

    %Mean wave energy dissipation
    delX = WGs(end)-WGs(1);
    ubr = sqrt(sum(ubj(end,:).^2));
    epsj = -1*(F(end,:)-F(1,:))./delX;
    fej = (epsj)/((0.6)*1025*ubr.*ubj(end,:).^2);
    fer = sum(fej.*ubj(end,:).^2)/sum(ubj(end,:).^2);
    HR.fe(i,:) = fer;
    HR.F(i,:) = trapz(F,2);
end
% HR.fe(HR.fe>0.8)=0.15;

%% Load processed model data -- LR
LRmodelFree = dir([LRmodelDir '*freeSurf_V2*']);
LRmodelFree = LRmodelFree(~ismember({LRmodelFree.name},{'.','..'}));
LRmodelFree = {LRmodelFree.name};

LRmodelBathy = dir([LRbathyDir '*_bathy_0.stl']);
LRmodelBathy = LRmodelBathy(~ismember({LRmodelBathy.name},{'.','..'}));
LRmodelBathy = {LRmodelBathy.name};

modelNames = 'Models_allLR.csv';
fid = fopen([workingDir modelNames]);
header = textscan(fgetl(fid),repmat('%s',1,9),'delimiter',',');
data = textscan(fid,repmat('%f',1,9),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);clear data

disp('Loading LR models')
%Loop through model scenarios and load the data into a structure
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(LRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    FS = load([LRmodelDir LRmodelFree{modelID}]);
    fn = fieldnames(FS);
    h = modelInfo.h(i);
    LR.h(i) = h;
    LR.Tp(i) = modelInfo.Tp(i);
    F = zeros(12,64);
    ubj = zeros(12,64); %orbital velocity from linear wave theory

    %Get local h with bathy file
    whichBathy = find(contains(LRmodelBathy,sprintf('%dm',h)));
    bathy = stlread([LRbathyDir LRmodelBathy{whichBathy}]);
    [bathyX,sortID] = sort(bathy.Points(:,1));
    bathyY = bathy.Points(sortID,3);
    
    %Crop time series based on start/stop times in modelInfo
    times = FS.time;
    timeID = find(times >= modelInfo.tStart(i) & times <= modelInfo.tStop(i));
    
    %Update 09/30/21: Calculation of wave dissipation factor fe from
    %surface elevation spectra (Lowe et al. 2005)
    WGs = FS.waveGauges;
    WGorder = [1 6 7 8 9 10 11 12 13 2 3 4 5]; %wave gauges are out of numerical order in FS
    for j = 1:length(WGorder)
        eta = FS.(fn{WGorder(j)}).eta(timeID);
        eta = repmat(eta,1,15); %make the time series longer
        [SSeta,f] = pwelch(eta,256,[],[],8);
        [~,bathyID] = min(abs(bathyX-WGs(j)));
        localh = -1*bathyY(bathyID); %local h at wave gauge
        
        k = wavek(f,localh);
        kh = k*localh;
        c = sqrt(9.81.*tanh(kh))./sqrt(k);
        n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
        Cg = n.*c;
        cutoff = find(f>=2,1,'first');
        F(j,:) = 1025*9.81*SSeta(2:cutoff).*Cg(2:cutoff);
        
        %Estimate near-bottom orbital velocity with linear wave theory
        aj = sqrt(2*SSeta(2:cutoff));
        omegaj = 2*pi*f(2:cutoff);
        ubj(j,:) = (aj.*omegaj)./sinh(kh(2:cutoff));
    end
    %Mean wave energy dissipation
    delX = WGs(end)-WGs(1);
    ubr = sqrt(sum(ubj(end,:).^2));
    epsj = -1*(F(end,:)-F(1,:))./delX;
    fej = (epsj)/((0.6)*1025*ubr.*ubj(end,:).^2);
    fer = sum(fej.*ubj(end,:).^2)/sum(ubj(end,:).^2);
    LR.fe(i,:) = fer;
    LR.F(i,:) = trapz(F,2);
end
% LR.fe(LR.fe>0.5)=0.25;

ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);
xProfiles = 37:0.5:40;
cc = flipud(hex2rgb({'#032760';'#0057d2';'#57b0ff';'#73dcff'})); %A dark -> light blue colormap (based on Bath 112)
% cc = [32 80 255;134 217 255;255 196 0;213 0 0]./255; %A red orange light blue blue colormap (based on GMT Panopoly)
markers = {'^';'s';'o';'d'};
lineStyle = {'-';'--';'-.';':'};
xLin = linspace(0,4,10);
yLin = linspace(0,4,10);
plot(xLin,yLin,'--k','linewidth',1.5);hold on

sp(1) = subplot(121);
xs = repmat(linspace(-0.6,0.6,4)',1,5);
waterDepths = unique(modelInfo.h);
wavePeriods = unique(modelInfo.Tp);
for i = 1:length(waterDepths)
    modelID = find(HR.h == waterDepths(i));
    fNorm = abs(HR.fe(modelID));
%     fNorm = mean(HR.F(modelID,2:end)./HR.F(modelID,1),2);
    x = HR.Tp(modelID);
    y = fNorm;
        edges = linspace(0,24,6);
        bins = discretize(x,edges);
        [meanVal,maxVal,stDev] = grpstats(y,bins,{@mean, @max, @std});
        errorbar(wavePeriods'+xs(i,:),meanVal,stDev,...
            sprintf('%s%s',markers{i},lineStyle{i}),'color',cc(i,:),...
            'markerfacecolor',cc(i,:),'markersize',6,'linewidth',1.2,'capsize',0);
    hold on;
end

sp(2) = subplot(122);
xs = repmat(linspace(-0.6,0.6,4)',1,5);
waterDepths = unique(modelInfo.h);
wavePeriods = unique(modelInfo.Tp);
for i = 1:length(waterDepths)
    modelID = find(LR.h == waterDepths(i));
    fNorm = abs(LR.fe(modelID));
%     fNorm = mean(LR.F(modelID,2:end)./LR.F(modelID,1),2);
    x = LR.Tp(modelID);
    y = fNorm;
        edges = linspace(0,24,6);
        bins = discretize(x,edges);
        [meanVal,maxVal,stDev] = grpstats(y,bins,{@mean, @max, @std});
        errorbar(wavePeriods'+xs(i,:),meanVal,stDev,...
            sprintf('%s%s',markers{i},lineStyle{i}),'color',cc(i,:),...
            'markerfacecolor',cc(i,:),'markersize',6,'linewidth',1.2,'capsize',0);
    hold on;
end

hModel = unique(HR.h);
%Make some dummy data for the legends
xs = repmat(linspace(1E5,1E6,2),length(hModel),1);
ys = repmat(linspace(min(x),max(x),2),length(hModel),1);
pp = zeros(length(hModel),1);
legText = cell(length(hModel),1);
for i = 1:length(hModel)
    pp(i) = plot(xs(i,:),ys(i,:),sprintf('%s%s',lineStyle{i},markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',6,'linewidth',1.2);hold on;
    legText{i} = sprintf('$h = %0.0f$ m',hModel(i));
end
leg = legend(pp,legText,'location','southeast','box','off','interpreter','latex');

%Global adjustments
set(sp(1),'xtick',4:4:20,'xlim',[2 22],...
    'ylim',[0 1.2],'ytick',0.2:0.2:1)
set(sp(2),'xtick',4:4:20,'xlim',[2 22],...
    'ylim',[0 1.2],'ytick',0:0.2:1,...
    'yticklabel',[])

%Positoning
set(sp(1),'position',[0.1 0.12 0.38 0.8])
set(sp(2),'position',[0.52 0.12 0.38 0.8])
set(leg,'position',[0.78 0.75 0.1 0.1])

%Labeling
ylabel(sp(1),'$\langle F\rangle / F_0 \   \mathrm{(-)}$','interpreter','latex')
xlabel(sp(1),'$T_p \ \mathrm{(s)}$','interpreter','latex')
xlabel(sp(2),'$T_p \ \mathrm{(s)}$','interpreter','latex')
title(sp(1),'High Relief')
title(sp(2),'Low Relief')
 
prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.008 0.008])
 
set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% print(ff,[workingDir 'Models_allHR_LR_feVsTp_alt1'],'-dpdf','-r0')
