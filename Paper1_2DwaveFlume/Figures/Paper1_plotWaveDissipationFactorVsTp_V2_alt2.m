% Create a plot of wave dissipation factor fe versus wave period for
% four different water depths
%
% Update:
% 10/8/21: Based on results in Gon et al. 2020, decided to plot fe vs.
% Ab/Kw (wave orbital excursion vs roughness) to normalize results. This
% figure will go in the discussion section of the paper. 
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
    omegaj = zeros(12,64);
    
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
        omegaj(j,:) = 2*pi*f(2:cutoff)';
        ubj(j,:) = (aj.*omegaj(j,:)')./sinh(kh(2:cutoff));  
    end

    %Mean wave energy dissipation
    delX = WGs(end)-WGs(1);
    ubr = sqrt(sum(ubj(end,:).^2));
    omegar = sum(omegaj(end,:).*(ubj(end,:).^2))/sum((ubj(end,:).^2));
    HR.Ab(i,:) = ubr/omegar;
    epsj = -1*(F(end,:)-F(1,:))./delX;
    fej = 4*(epsj)/(1025*ubr.*ubj(end,:).^2);
    fer = sum(fej.*ubj(end,:).^2)/sum(ubj(end,:).^2);
    HR.fe(i,:) = fer;
    HR.F(i,:) = trapz(F,2);
end
HR.Ab(HR.fe<0.005) = [];
HR.Tp(HR.fe<0.005) = [];
HR.fe(HR.fe<0.005) = [];

% Load processed model data -- LR
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
    omegaj = zeros(12,64);
    
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
        omegaj(j,:) = 2*pi*f(2:cutoff)';
        ubj(j,:) = (aj.*omegaj(j,:)')./sinh(kh(2:cutoff));
    end
    %Mean wave energy dissipation
    delX = WGs(end)-WGs(1);
    ubr = sqrt(sum(ubj(end,:).^2));
    omegar = sum(omegaj(end,:).*(ubj(end,:).^2))/sum((ubj(end,:).^2));
    LR.Ab(i,:) = ubr/omegar;
    epsj = -1*(F(end,:)-F(1,:))./delX;
    fej = 4*(epsj)/(1025*ubr.*ubj(end,:).^2);
    fer = sum(fej.*ubj(end,:).^2)/sum(ubj(end,:).^2);
    LR.fe(i,:) = fer;
    LR.F(i,:) = trapz(F,2);
end
LR.Ab(LR.fe<0.005) = [];
LR.Tp(LR.fe<0.005) = [];
LR.fe(LR.fe<0.005) = [];

ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   1000   400]);
xProfiles = 37:0.5:40;
cc = hex2rgb({'#1e5cb3';'#04a1e6';'#ffdf28';'#f58b35';'#bf171f'}); %blue-yellow-red colormap based on GMT panopoly
% cc = [32 80 255;134 217 255;255 196 0;213 0 0]./255; %A red orange light blue blue colormap (based on GMT Panopoly)
lineStyle = {'-';'--';'-.';':';'-'};
markers = {'o';'s';'d';'^';'p'};

sp(1) = subplot(131);
%Fit Nielsen (1992)
kw = 0.36; %from Paper1_calculateRugosity.m
Ab = [0.1 30]; %range of Ab from data
fe = exp(5.5.*((kw./Ab).^0.2)-6.3);
% fe = 4.2*((Ab./kw).^-1.57);
pfit1(2) = plot(Ab./kw,fe,'--k','linewidth',1);hold on

x = LR.Ab./kw;
y = abs(LR.fe);
B = polyfit(log(x),log(y),1);
Yfit = polyval(B,log(x));
Yfit2 = exp(Yfit);
pfit2(2) = plot(x,Yfit2,'-k','linewidth',1.2);
rsq = corrcoef(log(x),log(y));rsq = abs(rsq(1,2));fprintf('R^2 of LR model power fit is %0.2f\n',rsq);
lineEqn = sprintf('$f_e = %0.2f(A_b/k_w)^{%0.2f}$',B(2),B(1));

wavePeriods = unique(modelInfo.Tp);
pp = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    modelID = find(LR.Tp == wavePeriods(i));
    x = LR.Ab(modelID)./kw;
    y = abs(LR.fe(modelID));
    pp(i) = plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end

%Make best fit line legends
leg1 = legend([pp; pfit2(2)],[legText; lineEqn],'location','southwest','box','off','interpreter','latex');

sp(2) = subplot(132);
%Fit Nielsen (1992)
kw = 0.58; %from Paper1_calculateRugosity.m
Ab = [0.1 40]; %range of Ab from data
fe = exp(5.5.*((kw./Ab).^0.2)-6.3);
pfit1(1) = plot(Ab./kw,fe,'--k','linewidth',1);hold on

x = HR.Ab./kw;
y = abs(HR.fe);
B = polyfit(log(x),log(y),1);
Yfit = polyval(B,log(x));
Yfit2 = exp(Yfit);
pfit2(1) = plot(x,Yfit2,'-k','linewidth',1.2);
rsq = corrcoef(log(x),log(y));rsq = abs(rsq(1,2));fprintf('R^2 of HR model power fit is %0.2f\n',rsq);
lineEqn = sprintf('$f_e = %0.2f(A_b/k_w)^{%0.2f}$',B(2),B(1));

wavePeriods = unique(modelInfo.Tp);
pp = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    modelID = find(HR.Tp == wavePeriods(i));
    x = HR.Ab(modelID)./kw;
    y = abs(HR.fe(modelID));
    pp(i) = plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end

%Make best fit line legends
leg2 = legend([pp; pfit2(1)],[legText; lineEqn],'location','southwest','box','off','interpreter','latex');

sp(3) = subplot(133);
kw1 = 0.36; %from Paper1_calculateRugosity.m
kw2 = 0.58; %from Paper1_calculateRugosity.m
wavePeriods = unique(modelInfo.Tp);
pp = zeros(length(wavePeriods),1);
legText = cell(length(wavePeriods),1);
for i = 1:length(wavePeriods)
    hold on
    if i == 1
        plot([-50 50],[0 0],'k','linewidth',1)
        plot([0 0],[-10 10],'k','linewidth',1)
    end
    modelID = find(HR.Tp == wavePeriods(i));
    x = (LR.Ab(modelID)./kw1)-(HR.Ab(modelID)./kw2);
    y = abs(LR.fe(modelID))-abs(HR.fe(modelID));
    pp(i) = plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);
    legText{i} = sprintf('$T_p = %0.0f$ s',wavePeriods(i));
end

leg3 = legend(pp,legText,'location','southwest','box','off','interpreter','latex');
% 
% 
%Global adjustments
set(sp(1),'xlim',[0.35 50],...
    'ylim',[0.001 10])
set(sp(2),'xlim',[0.35 50],...
    'ylim',[0.001 10],...
    'yticklabel',[])
set([sp(1) sp(2)],'xscale','log','yscale','log')
set(sp(3),'xlim',[-22 40],...
    'ylim',[-4 4])

%Positoning
set(sp(1),'position',[0.08 0.12 0.27 0.8])
set(sp(2),'position',[0.37 0.12 0.27 0.8],'yticklabel',[])
set(sp(3),'position',[0.7 0.12 0.27 0.8])
set(leg1,'position',[0.105 0.185 0.15 0.15])
set(leg2,'position',[0.39 0.185 0.15 0.15])
set(leg3,'position',[0.665 0.195 0.15 0.15])

%Labeling
ylabel(sp(1),'$f_e \ \mathrm{(-)}$','interpreter','latex')
ylabel(sp(3),'${f_e}_{LR} - {f_e}_{HR} \ \mathrm{(-)}$','interpreter','latex')
xlabel(sp(1),'$A_b/k_w \ \mathrm{(-)}$','interpreter','latex')
xlabel(sp(2),'$A_b/k_w \ \ \mathrm{(-)}$','interpreter','latex')
xlabel(sp(3),'${A_b/k_w}_{LR} - {A_b/k_w}_{HR} \ \ \mathrm{(-)}$','interpreter','latex')
title(sp(1),'Low Relief')
title(sp(2),'High Relief')
title(sp(3),'Difference')
 
prettyfigures('text',11,'labels',12,'box',1,'tickdir','in','tlength',[0.008 0.008])
 
set(ff,'units','inches','renderer','painters');
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(ff,[workingDir 'Models_allHR_LR_feVsAbKw_V2_alt2'],'-dpdf','-r0')
