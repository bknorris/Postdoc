% Create reduced-size data files to simplify figure making. This script
% loads modeled data from the \DataAnalysis folder for HR and LR models and
% computes key variables (wave energy dissipation between the model inlet
% and patch, D; total mean TKE dissipation rate (eps); wave breaking
% parameter (lambda); wave steepness (steepness)). TKE dissipation rate is
% binned (averaged) along the x-z axes into 100 bins to reduce storage
% space.
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
versionNO = 'allHR_LR_V6';

%Define data paths
workingDir = 'c:\Users\bknorris\Documents\Models\Paper2_OptimizingRestoration\ModelRuns\Scenarios\';
HRmodelDir = 'c:\Users\bknorris\Documents\Models\Paper2_OptimizingRestoration\ModelRuns\Scenarios\DataAnalysis\';
% LRmodelDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\LR_Domain\';
% HRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\Models\';
% LRbathyDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\LR_Domain\Models\';

%Get modelNames from the model directory
HRmodelRaw = dir([HRmodelDir '*rawData_V1*']);
HRmodelRaw = HRmodelRaw(~ismember({HRmodelRaw.name},{'.','..'}));
HRmodelRaw = {HRmodelRaw.name};

HRmodelFree = dir([HRmodelDir '*freeSurf_V1*']);
HRmodelFree = HRmodelFree(~ismember({HRmodelFree.name},{'.','..'}));
HRmodelFree = {HRmodelFree.name};

% HRmodelBathy = dir([HRbathyDir '*_bathy_0.stl']);
% HRmodelBathy = HRmodelBathy(~ismember({HRmodelBathy.name},{'.','..'}));
% HRmodelBathy = {HRmodelBathy.name};
% 
% LRmodelRaw = dir([LRmodelDir '*rawData_V2*']);
% LRmodelRaw = LRmodelRaw(~ismember({LRmodelRaw.name},{'.','..'}));
% LRmodelRaw = {LRmodelRaw.name};
% 
% LRmodelFree = dir([LRmodelDir '*freeSurf_V2*']);
% LRmodelFree = LRmodelFree(~ismember({LRmodelFree.name},{'.','..'}));
% LRmodelFree = {LRmodelFree.name};
% 
% LRmodelBathy = dir([LRbathyDir '*_bathy_0.stl']);
% LRmodelBathy = LRmodelBathy(~ismember({LRmodelBathy.name},{'.','..'}));
% LRmodelBathy = {LRmodelBathy.name};

disp('Loading HR models')
t1 = tic;
%Load modelInfo from the CSV file and create a structure from the data
modelNames = 'adjustModelSampling_TEST.csv';
fid = fopen([workingDir modelNames]);
header = textscan(fgetl(fid),repmat('%s',1,8),'delimiter',',');
data = textscan(fid,repmat('%f',1,8),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);clear data
%Loop through model scenarios and load the data into a structure
data = struct();
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.orgScenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(HRmodelRaw','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.orgScenarioNumber),HRmodelRaw{modelID})
    HR = load([HRmodelDir HRmodelRaw{modelID}]);
%     FS = load([HRmodelDir HRmodelFree{modelID}]);
    
    Tp = modelInfo.wavePeriod(i);
    Hs = modelInfo.waveHeight(i);
    h = modelInfo.waterDepth(i);
    gamma = Hs/h;
    lambda = (9.81.*(Tp.^2))./(2*pi);
    steepness = Hs/lambda;
%     fn = fieldnames(FS);
%     F = zeros(12,64);
%     Sub = zeros(12,64); %orbital velocity from linear wave theory
%     ub = zeros(12,1);
%     
%     %Get local h with bathy file
%     whichBathy = find(contains(HRmodelBathy,sprintf('%dm',h)));
%     bathy = stlread([HRbathyDir HRmodelBathy{whichBathy}]);
%     [bathyX,sortID] = sort(bathy.Points(:,1));
%     bathyY = bathy.Points(sortID,3);
%     
%     %% Update: 10/4/21: Calculation of wave energy flux, dissipation, and
%     %energy dissipation factor (fe)
%     WGs = FS.waveGauges;
%     WGorder = [1 6 7 8 9 10 11 12 13 2 3 4 5]; %wave gauges are out of numerical order in FS
%     for j = 1:length(WGorder)
%         [SSeta,f] = pwelch(FS.(fn{WGorder(j)}).eta,[],[],[],8);
%         [~,bathyID] = min(abs(bathyX-WGs(j)));
%         localh = -1*bathyY(bathyID); %local h at wave gauge
%         
%         k = wavek(f,localh);
%         kh = k*localh;
%         c = sqrt(9.81.*tanh(kh))./sqrt(k);
%         n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
%         Cg = n.*c;
%         cutoff = find(f>=2,1,'first');
%         F(j,:) = 1025*9.81*SSeta(2:cutoff).*Cg(2:cutoff);
%         
%         %Estimate near-bottom orbital velocity with linear wave theory
%         %(e.g., Wiberg & Sherwood 2008, eq. 5-6)
%         Kub = (2*pi*f(2:cutoff))./sinh(kh(2:cutoff));
%         Sub(j,:) = (abs(Kub).^2).*SSeta(2:cutoff); 
%         ub(j) = sqrt(2*trapz(Sub(j,:))*mean(diff(f))); 
%     end
%     %Mean wave energy dissipation
%     delX = WGs(end)-WGs(2);
%     delF = F(end,:)-F(2,:);
%     D = -1*(delF./delX);
%     fe = (2*D)\(1025*(sqrt(2/pi))*ub(end).*Sub(end,:));
        
    %% Average turbulence in time 
    %use tStart and tStop to crop model data
    %based on analysis of "good" timesteps (see
    %Paper1_animateModelTimesteps.m)
%     if min(size(HR.epsilon.epsilon)) ~= length(HR.time) %some cases have epsilon shorter than U and time variables, catch these
%         times = HR.time(1:min(size(HR.epsilon.epsilon)));
%     else
%         times = HR.time;
%     end
%     timeID = find(times >= modelInfo.tStart(i) & times <= modelInfo.tStop(i));
    epsAvg = mean(HR.epsilon.epsilon(:,:),2);
    x = HR.epsilon.x(:,1);
    z = HR.epsilon.z(:,1);
    
    nBins = 100;
    xBins = linspace(min(x),max(x),100);
    zBins = linspace(min(z),max(z),34);
%     eBins = linspace(min(epsAvg),max(epsAvg),nBins);
    
    ix = discretize(x, xBins);
    iz = discretize(z, zBins);
    idx = sub2ind([100 34], ix, iz);
    epsg = accumarray(idx, epsAvg, [100*34 1], @(x) mean(x), NaN);
    epsg = reshape(epsg, 100, 34);
    
    %Compute Umag and average in time
    Umag = sqrt((HR.U.Ux.^2)+(HR.U.Uy.^2)+(HR.U.Uz.^2));
    UmagAvg = mean(Umag(:,:),2);
    Umagg = accumarray(idx, UmagAvg, [nBins*nBins 1], @(x) mean(x), NaN);
    Umagg = reshape(Umagg, nBins, nBins);
    
    %Average TKE (k) in time
    TKEavg = mean(HR.k.k(:,:),2);
    TKEavgg = accumarray(idx, TKEavg, [nBins*nBins 1], @(x) mean(x), NaN);
    TKEavgg = reshape(TKEavgg, nBins, nBins);
    
    %% Compute mean oscillatory velocity Uw and orbital velocity Ubr
    %Uw^3 should be proportional to TKE dissipation from bottom 
    %friction if dissipation comes from bottom friction and not breaking
    [m,n] = size(HR.U.x(:,timeID));
    Ux = zeros(nBins,nBins,n);
    Uy = zeros(nBins,nBins,n);
    Uz = zeros(nBins,nBins,n);
    for j = 1:n
        int = accumarray(idx, HR.U.Ux(:,timeID(j)), [nBins*nBins 1], @(x) mean(x), NaN);
        Ux(:,:,j) = reshape(int, nBins, nBins);
        int = accumarray(idx, HR.U.Uy(:,timeID(j)), [nBins*nBins 1], @(x) mean(x), NaN);
        Uy(:,:,j) = reshape(int, nBins, nBins);
        int = accumarray(idx, HR.U.Uz(:,timeID(j)), [nBins*nBins 1], @(x) mean(x), NaN);
        Uz(:,:,j) = reshape(int, nBins, nBins);
    end
    %Calculate RMS velocities (Luhar et al. 2013)
    %Compute wave orbital velocities (Ubr) for all data points.
    Ubr = zeros(nBins,nBins);
    Uc = zeros(nBins,nBins);
    Uw = zeros(nBins,nBins);
    Ww = zeros(nBins,nBins);
    for j = 1:nBins
        for k = 1:nBins
            Uxx =  squeeze(Ux(j,k,:));
            Uyy =  squeeze(Uy(j,k,:));
            Uzz =  squeeze(Uz(j,k,:));
            if any(isnan(Uxx) | isnan(Uyy) | isnan(Uzz))
                Ubr(j,k) = NaN;
                Uc(j,k) = NaN;
                Uw(j,k) = NaN;
                Ww(j,k) = NaN;
            else
                [Cuu,f] = pwelch(Uxx,[],[],[],8);
                ff = find(f>=0.03 & f<=1);
                Ubr(j,k) = sqrt(2*trapz(Cuu(ff)*mean(diff(f))));
                Ec = (1/n)*sum(Uxx);Nc = (1/n)*sum(Uyy);
                Ewrms = sqrt((1/n)*sum((Uxx-Ec).^2));
                Nwrms = sqrt((1/n)*sum((Uyy-Nc).^2));
                Wc = (1/n)*sum(Uzz);Wwrms = sqrt((1/n)*sum((Uzz-Wc).^2));
                Uc(j,k) = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
                Uw(j,k) = sqrt(2)*Uwrms;
                Ww(j,k) = sqrt(2)*Wwrms;
            end
        end
    end
    
    %Load data structure
    data.scenarioNumber(i) = str2double(whichModel);
    data.Tp(i) = Tp;
    data.Hs(i) = Hs;
    data.h(i) = h;
    data.gamma(i) = gamma;
    data.steepness(i) = steepness;
%     data.F(:,i) = F;
%     data.D(i) = D;
%     data.fe(i) = fe;
    data.xBins(:,:,i) = xBins;
    data.zBins(:,:,i) = zBins;
    data.k(:,:,i) = TKEavgg;
    data.eps(:,:,i) = epsg;
    data.Umag(:,:,i) = Umagg;
    data.Ubr(:,:,i) = Ubr;
    data.ub(:,i) = ub;
    data.Uc(:,:,i) = Uc;
    data.Uw(:,:,i) = Uw;
    data.Ww(:,:,i) = Ww;
end
save([HRmodelDir 'HR_postProcess_epsilon_' versionNO '.mat'],'-struct','data','-mat','-v7.3')
clear data
disp(['Data analysis completed in ' num2str(toc(t1)/60) ' minutes'])

disp('Loading LR models')
t2 = tic;
%Load modelInfo from the CSV file and create a structure from the data
modelNames = 'Models_allLR.csv';
fid = fopen([workingDir modelNames]);
header = textscan(fgetl(fid),repmat('%s',1,9),'delimiter',',');
data = textscan(fid,repmat('%f',1,9),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);clear data
%Loop through model scenarios and load the data into a structure
data = struct();
for i = 1:length(modelInfo.scenarioNumber)
    whichModel = string(modelInfo.scenarioNumber(i));
    %Find model folders
    scenarioNumber = regexp(LRmodelFree','(?<=_)(.*?)(?=_)','match');
    scenarioNumber = cat(1, scenarioNumber{:});
    modelID = find(strcmp(whichModel,scenarioNumber));
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelInfo.scenarioNumber),LRmodelFree{modelID})
    LR = load([LRmodelDir LRmodelRaw{modelID}]);
    FS = load([LRmodelDir LRmodelFree{modelID}]);
    
    Tp = modelInfo.Tp(i);
    Hs = modelInfo.Hs(i);
    h = modelInfo.h(i);
    gamma = Hs/h;
    lambda = (9.81.*(Tp.^2))./(2*pi);
    steepness = Hs/lambda;
    fn = fieldnames(FS);
    F = zeros(12,1);
    ub = zeros(12,1); %orbital velocity from linear wave theory
    
    %Get local h with bathy file
    whichBathy = find(contains(LRmodelBathy,sprintf('%dm',h)));
    bathy = stlread([LRbathyDir LRmodelBathy{whichBathy}]);
    [bathyX,sortID] = sort(bathy.Points(:,1));
    bathyY = bathy.Points(sortID,3);
    
    %% Update: 10/4/21: Calculation of wave energy flux, dissipation, and
    %energy dissipation factor (fe)
    WGs = FS.waveGauges;
    WGorder = [1 6 7 8 9 10 11 12 13 2 3 4 5]; %wave gauges are out of numerical order in FS
    for j = 1:length(WGorder)-1
        [SSeta,f] = pwelch(FS.(fn{WGorder(j)}).eta,[],[],[],8);
        [~,bathyID] = min(abs(bathyX-WGs(j)));
        localh = -1*bathyY(bathyID); %local h at wave gauge
        
        k = wavek(f,localh);
        kh = k*localh;
        c = sqrt(9.81.*tanh(kh))./sqrt(k);
        n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
        Cg = n.*c;
        cutoff = find(f>=2,1,'first');
        F(j) = trapz(1025*9.81*SSeta(2:cutoff).*Cg(2:cutoff));
        
        %Estimate near-bottom orbital velocity with linear wave theory
        %(e.g., Wiberg & Sherwood 2008, eq. 5-6)
        Kub = (2*pi*f(2:cutoff))./sinh(kh(2:cutoff));
        Sub = (abs(Kub).^2).*SSeta(2:cutoff); 
        ub(j) = sqrt(2*trapz(Sub)*mean(diff(f))); 
    end
    %Mean wave energy dissipation
    delX = WGs(end)-WGs(2);
    delF = F(end)-F(2);
    D = -1*(delF/delX);
    fe = D/((2/(3*pi))*1025*(mean(ub(2:end))^3));
        
    %% Average turbulence in time 
    %use tStart and tStop to crop model data
    %based on analysis of "good" timesteps (see
    %Paper1_animateModelTimesteps.m)
    if min(size(LR.epsilon.epsilon)) ~= length(LR.time) %some cases have epsilon shorter than U and time variables, catch these
        times = LR.time(1:min(size(LR.epsilon.epsilon)));
    else
        times = LR.time;
    end
    timeID = find(times >= modelInfo.tStart(i) & times <= modelInfo.tStop(i));
    epsAvg = mean(LR.epsilon.epsilon(:,timeID),2);
    x = LR.epsilon.x(:,1)-30;
    z = LR.epsilon.z(:,1);
    
    nBins = 100;
    xBins = linspace(min(x),max(x),nBins);
    zBins = linspace(min(z),max(z),nBins);
    eBins = linspace(min(epsAvg),max(epsAvg),nBins);
    
    ix = discretize(x, xBins);
    iz = discretize(z, zBins);
    idx = sub2ind([nBins nBins], ix, iz);
    epsg = accumarray(idx, epsAvg, [nBins*nBins 1], @(x) mean(x), NaN);
    epsg = reshape(epsg, nBins, nBins);
    
    %Compute Umag and average in time
    Umag = sqrt((LR.U.Ux.^2)+(LR.U.Uy.^2)+(LR.U.Uz.^2));
    UmagAvg = mean(Umag(:,timeID),2);
    Umagg = accumarray(idx, UmagAvg, [nBins*nBins 1], @(x) mean(x), NaN);
    Umagg = reshape(Umagg, nBins, nBins);
    
    %Average TKE (k) in time
    TKEavg = mean(LR.k.k(:,timeID),2);
    TKEavgg = accumarray(idx, TKEavg, [nBins*nBins 1], @(x) mean(x), NaN);
    TKEavgg = reshape(TKEavgg, nBins, nBins);
    
    %% Compute mean oscillatory velocity Uw and orbital velocity Ubr
    %Uw^3 should be proportional to TKE dissipation from bottom 
    %friction if dissipation comes from bottom friction and not breaking
    [m,n] = size(LR.U.x(:,timeID));
    Ux = zeros(nBins,nBins,n);
    Uy = zeros(nBins,nBins,n);
    Uz = zeros(nBins,nBins,n);
    for j = 1:n
        int = accumarray(idx, LR.U.Ux(:,timeID(j)), [nBins*nBins 1], @(x) mean(x), NaN);
        Ux(:,:,j) = reshape(int, nBins, nBins);
        int = accumarray(idx, LR.U.Uy(:,timeID(j)), [nBins*nBins 1], @(x) mean(x), NaN);
        Uy(:,:,j) = reshape(int, nBins, nBins);
        int = accumarray(idx, LR.U.Uz(:,timeID(j)), [nBins*nBins 1], @(x) mean(x), NaN);
        Uz(:,:,j) = reshape(int, nBins, nBins);
    end
    %Calculate RMS velocities (Luhar et al. 2013)
    %Compute wave orbital velocities (Ubr) for all data points.
    Ubr = zeros(nBins,nBins);
    Uc = zeros(nBins,nBins);
    Uw = zeros(nBins,nBins);
    Ww = zeros(nBins,nBins);
    for j = 1:nBins
        for k = 1:nBins
            Uxx =  squeeze(Ux(j,k,:));
            Uyy =  squeeze(Uy(j,k,:));
            Uzz =  squeeze(Uz(j,k,:));
            if any(isnan(Uxx) | isnan(Uyy) | isnan(Uzz))
                Ubr(j,k) = NaN;
                Uc(j,k) = NaN;
                Uw(j,k) = NaN;
                Ww(j,k) = NaN;
            else
                [Cuu,f] = pwelch(Uxx,[],[],[],8);
                ff = find(f>=0.03 & f<=1);
                Ubr(j,k) = sqrt(2*trapz(Cuu(ff)*mean(diff(f))));
                Ec = (1/n)*sum(Uxx);Nc = (1/n)*sum(Uyy);
                Ewrms = sqrt((1/n)*sum((Uxx-Ec).^2));
                Nwrms = sqrt((1/n)*sum((Uyy-Nc).^2));
                Wc = (1/n)*sum(Uzz);Wwrms = sqrt((1/n)*sum((Uzz-Wc).^2));
                Uc(j,k) = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
                Uw(j,k) = sqrt(2)*Uwrms;
                Ww(j,k) = sqrt(2)*Wwrms;
            end
        end
    end
    
    %Load data structure
    data.scenarioNumber(i) = str2double(whichModel);
    data.Tp(i) = Tp;
    data.Hs(i) = Hs;
    data.h(i) = h;
    data.gamma(i) = gamma;
    data.steepness(i) = steepness;
    data.F(:,i) = F;
    data.D(i) = D;
    data.fe(i) = fe;
    data.xBins(:,:,i) = xBins;
    data.zBins(:,:,i) = zBins;
    data.k(:,:,i) = TKEavgg;
    data.eps(:,:,i) = epsg;
    data.Umag(:,:,i) = Umagg;
    data.Ubr(:,:,i) = Ubr;
    data.ub(:,i) = ub;
    data.Uc(:,:,i) = Uc;
    data.Uw(:,:,i) = Uw;
    data.Ww(:,:,i) = Ww;
end
save([LRmodelDir 'LR_postProcess_epsilon_' versionNO '.mat'],'-struct','data','-mat','-v7.3')
clear data
disp(['Data analysis completed in ' num2str(toc(t2)/60) ' minutes'])

