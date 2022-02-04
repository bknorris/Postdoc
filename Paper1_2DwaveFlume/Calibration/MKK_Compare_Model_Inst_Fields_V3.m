%Compare model derived quantities (e.g., U, epsilon) with AQDP measurements.
%
% Updates:
% 01/13/21: Use raw data files instead of wavestats for flexibility, plot
% burst-mean oscillatory velocities instead of mean velocities
% 01/29/21: Revised script to plot results from multiple bursts together
%
% This is version 3 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
lineName = {'line2_alpha.water_k_nut_omega_p_rgh';'line2_U'};
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
insts = 'MKK18HR101aqdHR-b.nc';
modelName = 'MKK_CombinedHR_V5';
% RBRbursts = [151 158]; %RBR bursts to analyze (range: 141 - 162)
% AQDbursts = [61 68]; %ADCP bursts to analyze (range: 51 - 72)
RBRbursts = [160]; %RBR bursts to analyze (range: 141 - 162)
AQDbursts = [70]; %ADCP bursts to analyze (range: 51 - 72)
saveFigs = 0; %save figures to disk; 0 is off, 1 is on

%% Load and process the model data
if length(RBRbursts) == 1
    burstTxt = sprintf('%d',RBRbursts(1));
else
    burstTxt = sprintf('%d-%d',RBRbursts(1),RBRbursts(2));
end
modelBasePath = ['d:\USGS\Models\MKK_Combined\2D\Burst' burstTxt '\' modelName '\MKKmodel\postProcessing\line\'];
% modelBasePath = ['c:\Users\user\Documents\Models\MKK_Combined\2D\Burst' burstTxt '\' modelName '\MKKmodel\postProcessing\line\'];
modelFolders = dir(modelBasePath);
modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));

profileLength = 81;
Em = zeros(profileLength,length(modelFolders)-2); %model epsilon
Um = zeros(profileLength,length(modelFolders)-2); %model U
Vm = zeros(profileLength,length(modelFolders)-2); %model V
Wm = zeros(profileLength,length(modelFolders)-2); %model W
for i = 2:2684 %length(modelFolders)-1
    fileName = dir([modelBasePath modelFolders(i).name '\']);
    fileName = {fileName.name};
    %Load epsilon data
    isFile = contains(fileName,lineName{1});
    fid = fopen([modelBasePath modelFolders(i).name '\' fileName{isFile}]);
    if contains(fileName{isFile},'epsilon')
        data = textscan(fid,'%n%n%n%n');
        freeS = find(data{2} >= 0.5);
        %Determine if file is epsilon or omega (closure scheme model)
        Em(1:length(freeS),i-1) = data{3}(freeS);
    else
        data = textscan(fid,'%n%n%n%n%n%n');
        freeS = find(data{2} >= 0.5);
        %Determine if file is epsilon or omega (closure scheme model)
        k = data{3}(freeS);
        w = data{5}(freeS);
        Cmu = 0.09; %Constant
        Em(1:length(freeS),i-1) = Cmu.*k.*w;
    end
    fclose(fid);
    clear data
    
    %Load U data
    isFile = contains(fileName,lineName{2});
    fid = fopen([modelBasePath modelFolders(i).name '\' fileName{isFile}]);
    data = textscan(fid,'%n%n%n%n');
    Um(1:length(freeS),i-1) = data{2}(freeS);
    Vm(1:length(freeS),i-1) = data{3}(freeS);
    Wm(1:length(freeS),i-1) = data{4}(freeS);
    if i == 3
        depthM = data{1};
    end
    fclose(fid);
    clear data
end

%Compute wave statistics from model data
Uwrms_m = zeros(profileLength,1); %model wave velocities
Uw_m = zeros(profileLength,1);
Wwrms_m = zeros(profileLength,1);
Ww_m = zeros(profileLength,1);
E_m = zeros(profileLength,1);
for i = 1:profileLength
    U = Um(i,1:2684);
    V = Vm(i,1:2684);
    W = Wm(i,1:2684);
    nwin = length(U);
    
    %Calculate RMS velocities (Lowe et al. 2005; Luhar et al. 2013)
    %Ec and Nc are mean east and north current velocities
    %Ewrms and Nwrms are rms oscillatory velocities
    %Uw and Ww are total mean oscillatory horizontal and vertical velocities
    Ec = (1/nwin)*sum(U);Nc = (1/nwin)*sum(V);Wc = (1/nwin)*sum(W);
    Ewrms = sqrt((1/nwin)*sum((U-Ec).^2));
    Nwrms = sqrt((1/nwin)*sum((V-Nc).^2));
    Wwrms_m(i) = sqrt((1/nwin)*sum((W-Wc).^2));
    Uwrms_m(i) = sqrt(Ewrms^2+Nwrms^2);
    Uw_m(i) = sqrt(2)*Uwrms_m(i);
    Ww_m(i) = sqrt(2)*Wwrms_m(i);
    E_m(i) = mean(Em(i,1:2684));
end

%% Load and process the instrument data
data = loadnc([dataPath insts]);

%Extract metadata from the netCDF file
fs = data.Gatts.instmeta_AQDSamplingrate; %sample rate
zp = data.Gatts.transducer_offset_from_bottom; %height of pressure sensor
zuv = zp - data.bindist; %height of velocity bin above bed
cellSize = data.Gatts.instmeta_AQDCellSize/1000; %cell size in mm

%Define working variables
if length(AQDbursts) == 1
    idx = find(data.burst == AQDbursts(1));
else
    idx = find(data.burst >= AQDbursts(1) & data.burst <= AQDbursts(2));
end
[m,n] = size(data.dn(:,idx));
P = reshape(data.P_1(:,idx),m*n,1)+zp;
Ui = reshape(data.u_1205(:,:,idx)./100,19,m*n); %convert cm/s to m/s;
Vi = reshape(data.v_1206(:,:,idx)./100,19,m*n);
Wi = reshape(data.w_1204(:,:,idx)./100,19,m*n);

% Ui = mean(data.u_1205(:,:,idx)./100,3); %convert cm/s to m/s;
% Vi = mean(data.v_1206(:,:,idx)./100,3);
% Wi = mean(data.w_1204(:,:,idx)./100,3);

%Compute wave statistics per depth bin
Uwrms_i = zeros(length(data.depth),1); %observed wave velocities
Uw_i = zeros(length(data.depth),1);
Wwrms_i = zeros(length(data.depth),1);
Ww_i = zeros(length(data.depth),1);
for i = 1:length(data.depth)
    U = Ui(i,:);
    V = Ui(i,:);
    W = Ui(i,:);
    nwin = length(U);
    
    %Bridge short gaps in time series if NaNs are present
    if any([isnan(U) isnan(V) isnan(W)])
        disp('Found NaNs in Velocity... interpolating with cmgbridge')
        fprintf('Bin : %0.f\n',i)
        U = cmgbridge(U,100,100,1000);
        V = cmgbridge(V,100,100,1000);
        W = cmgbridge(W,100,100,1000);
    end
    
    %Calculate RMS velocities
    Ec = (1/nwin)*sum(U);Nc = (1/nwin)*sum(V);Wc = (1/nwin)*sum(W);
    Ewrms = sqrt((1/nwin)*sum((U-Ec).^2));
    Nwrms = sqrt((1/nwin)*sum((V-Nc).^2));
    Wwrms_i(i) = sqrt((1/nwin)*sum((W-Wc).^2));
    Uwrms_i(i) = sqrt(Ewrms^2+Nwrms^2);
    Uw_i(i) = sqrt(2)*Uwrms_i(i);
    Ww_i(i) = sqrt(2)*Wwrms_i(i);
end

%Compute dissipation rate of turbulence (Wiles et al., 2006)

r = (cellSize)/cosd(25);   %along-beam distance (AQDPs have beam angles of 25 deg)
lags = 5;
[m,n] = size(Wi);
maxbin = m-lags;
itt = 1;                   %iteration number
D = zeros(maxbin,lags);    %structure function D(z,r)
R = zeros(maxbin,lags);    %along-beam distance r^2/3
epsilon = zeros(1,m); %dissipation rate [m^2/s^3]
rr = zeros(1,lags);
while itt <= maxbin
    d = zeros(n,lags);
    idl = itt:itt+lags;
    for kk = 1:lags %compute velocity differences (vel(z,1) - vel(z,2:lags))^2
        d(:,kk) = (Wi(itt,:)-(Wi(idl(kk+1),:))).^2;
        D(itt,kk) = nanmean(d(:,kk));
        rr(:,kk) = r*kk;
    end
    R(itt,:) = rr.^(2/3);
    [pf,S] = polyfit(R(itt,:),D(itt,:),1);  %linear regression of D along r^2/3
    A = pf(1);
    epsilon(idl(4)) = (A/2.1)^(3/2);
    itt = itt+1;
end
epsilon(abs(imag(epsilon))>0) = NaN;

%% Find depths consistent between the model and the instrument
nearBed_m = find(depthM >= 0 & depthM <= max(zuv));

%Plot some simple comparisons for now
ff(1) = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[800 350   1250   500]);
sp(1) = subplot(131);
hold on
p(1) = plot(Uwrms_i,zuv,'linewidth',1.5,'color','b');
p(2) = plot(Uwrms_m(nearBed_m),depthM(nearBed_m),'linewidth',1.5,'color','r');
hold off
sp(2) = subplot(132);
hold on
p(1) = plot(Wwrms_i,zuv,'linewidth',1.5,'color','b');
p(2) = plot(Wwrms_m(nearBed_m),depthM(nearBed_m),'linewidth',1.5,'color','r');
hold off
sp(3) = subplot(133);
hold on
p(1) = plot(epsilon,zuv,'linewidth',1.5,'color','b');
p(2) = plot(E_m(nearBed_m),depthM(nearBed_m),'linewidth',1.5,'color','r');
hold off
leg = legend(p,{'AQDP';'Model'});
ylabel(sp(1),'\bf\itDepth [m]')
xlabel(sp(1),'\bf\itU_{w,RMS} [m/s]')
xlabel(sp(2),'\bf\itW_{w,RMS} [m/s]')
xlabel(sp(3),'\bf\it\epsilon [m^2/s^2]')
set(sp,'ylim',[0.05 0.7],'ytick',0:0.1:0.7)
set([sp(1) sp(2)],'xlim',[0 0.15])
set(sp(3),'xscale','log','xlim',[10^-7 10^-3])
set(sp(1),'position',[0.1 0.15 0.23 0.72])
set(sp(2),'position',[0.4 0.15 0.23 0.72])
set(sp(3),'position',[0.7 0.15 0.23 0.72])
set(leg,'position',[0.9 0.7 0.05 0.05])
suptitle(['\bf\it' sprintf('Burst no. %s',burstTxt)])
prettyfigures('text',11,'labels',13)
if saveFigs == 1
    export_fig(ff(1),[figPath modelName '_AQDP_MODEL-Fields_V3'],'-png','-nocrop')
end
