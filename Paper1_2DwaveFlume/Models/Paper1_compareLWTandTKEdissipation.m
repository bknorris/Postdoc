% Compare Wave energy dissipation from linear wave theory and energy
% dissipation estimated directly from the model
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
HRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\Models\';

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
    delX = WGs(end)-WGs(2);
    ubr = sqrt(sum(mean(ubj(2:end,:)).^2));
    omegar = sum(omegaj(end,:).*(mean(ubj(2:end,:)).^2))/sum((mean(ubj(2:end,:)).^2));
    HR.Ab(i,:) = ubr/omegar;
    epsj = -1*(F(end,:)-F(2,:))./delX;
    fej = 4*(epsj)/(1025*ubr.*mean(ubj(2:end,:)).^2);
    fer = sum(fej.*mean(ubj(2:end,:)).^2)/sum(mean(ubj(2:end,:)).^2);
    HR.ubr(i,:) = ubr;
    HR.eps(i,:) = sum(epsj);
end
FS = HR;
%Load processed model data
HR = load([HRmodelDir 'HR_postProcess_epsilon_allHR_LR_V4.mat']);


%Create figure
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);
lineStyle = {'-';'--';'-.';':';'-'};
markers = {'o';'s';'d';'^';'p'};
sp = zeros(2,1);
cc = hex2rgb({'#1e5cb3';'#04a1e6';'#ffdf28';'#f58b35';'#bf171f'}); %blue-yellow-red colormap based on GMT panopoly
% cc = flipud(hex2rgb({'#032760';'#0057d2';'#57b0ff';'#73dcff'})); %A dark -> light blue colormap (based on Bath 112)
%Do the HR models first; make sure these are the same as
%Paper1_plotDepthProfiles.m!
patchLo = 35;
patchHi = 41;
zPatch = [-0.9 -1.9 -2.9 -3.9];
upper = 15;
wavePeriods = unique(modelInfo.Tp);
for i = 1:length(wavePeriods)
    %Collect data by wave period then water depth!
    modelID = find(HR.Tp == wavePeriods(i));
    EPS = zeros(length(modelID),1);
    
    %EXPERIMENTAL: try filtering results by their proximity to the bed
    for j = 1:length(modelID)
        eps1 = squeeze(HR.eps(:,:,modelID(j)));
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
        data = eps1(xLo:xHi,:);[m,n] = size(data);
        uwfilt = NaN(m,n);
        for k = 1:m
            nanID = find(isnan(data(k,1:end-1)),1,'last');
            uwfilt(k,nanID+1:nanID+upper) = data(k,nanID+1:nanID+upper);
        end  
        EPS(j) = nanmean(nanmean(uwfilt));
    end
    modelID = find(FS.Tp == wavePeriods(i));
    eps = abs(FS.eps(modelID));
    x = eps;
    y = ((EPS*1025*0.4).^(2/3)).*(((3*pi)/8).*FS.ubr(modelID)); %boundary layer scaling from huang 2012
    plot(x,y,sprintf('%s',markers{i}),...
        'color',cc(i,:),'markerfacecolor',cc(i,:),'markersize',5,'linewidth',1);hold on
    (y./x)*100
end
