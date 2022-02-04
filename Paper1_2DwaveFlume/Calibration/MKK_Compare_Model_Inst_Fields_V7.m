%Compare model derived quantities (e.g., U, epsilon) with AQDP measurements.
%
% Updates:
% 01/13/21: Use raw data files instead of wavestats for flexibility, plot
% burst-mean oscillatory velocities instead of mean velocities
% 01/29/21: Revised script to plot results from multiple bursts together
% 02/08/21: Added functionality to plot multiple model results together on
% the same figure.
% 02/11/21: Added functionality to plot either multiple models OR multiple
% sample lines from the same model on the same figure.
% 03/10/21: Replaced multi line plot with averaged line plot (averages
% multiple lines together).
%
% This is version 7 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
modelNames = {'MKK_CombinedHR_V2'};
modelFigName = 'MKK_CombinedHR_V2';
insts = 'MKK18HR101aqdHR-b.nc';
instName = {'C15';'HR-Aquadopp'};
modelLegendText = {'Observations';'V27, Mesh [0.3 - 0.04 m]'};
lineName = {'line6_alpha.water';'line6_U'};

RBRbursts = [151 158]; %RBR bursts to analyze (range: 141 - 162)
AQDbursts = [61 68]; %ADCP bursts to analyze (range: 51 - 72)
% RBRbursts = [160]; %RBR bursts to analyze (range: 141 - 162)
% AQDbursts = [70]; %ADCP bursts to analyze (range: 51 - 72)
saveFigs = 0; %save figures to disk; 0 is off, 1 is on
saveData = 0; %save out profile data to disk; 0 is off, 1 is on
multipleModels = 0; %plot multiple models together; 0 is off, 1 is on
avgLines = 1; %average multiple lines together; 0 is off, 1 is on

%% Load and process the model data
profileLength = 1000;
if multipleModels == 1
    multiPlot = length(modelNames);
    ln = [1 2;1 2;1 2;3 4;1 2;5 6]; %specify which lines to use (must be length(modelNames))
    Em = NaN(profileLength,4096,multiPlot); %model epsilon
    Um = NaN(profileLength,4096,multiPlot); %model U
    Vm = NaN(profileLength,4096,multiPlot); %model V
    Wm = NaN(profileLength,4096,multiPlot); %model W
    depthM = NaN(profileLength,1,multiPlot);%model depth
    for i = 1:multiPlot
        if length(RBRbursts) == 1
            burstTxt = sprintf('%d',RBRbursts(1));
        else
            burstTxt = sprintf('%d-%d',RBRbursts(1),RBRbursts(2));
        end
        modelBasePath = ['c:\Users\user\Documents\Models\MKK_Combined\2D\Burst' burstTxt '\' modelNames{i} '\MKKmodel\postProcessing\line\'];
        modelFolders = dir(modelBasePath);
        modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
        for j = 2:length(modelFolders)-1
            fileName = dir([modelBasePath modelFolders(j).name '\']);
            fileName = {fileName.name};
            %Load epsilon data
            isFile = contains(fileName,lineName{ln(i,1)});
            fid = fopen([modelBasePath modelFolders(j).name '\' fileName{isFile}]);
            if contains(fileName{isFile},'epsilon')
                data = textscan(fid,'%n%n%n%n%n');
                freeS = find(data{2} >= 0.5);
                %Determine if file is epsilon or omega (closure scheme model)
                Em(1:length(freeS),j-1,i) = data{3}(freeS);
            else
                data = textscan(fid,'%n%n%n%n%n%n');
                freeS = find(data{2} >= 0.5);
                %Determine if file is epsilon or omega (closure scheme model)
                k = data{3}(freeS);
                w = data{4}(freeS);
                Cmu = 0.09; %Constant
                Em(1:length(freeS),j-1,i) = Cmu.*k.*w;
            end
            
            fclose(fid);
            clear data
            
            %Load U data
            isFile = contains(fileName,lineName{ln(i,2)});
            fid = fopen([modelBasePath modelFolders(j).name '\' fileName{isFile}]);
            data = textscan(fid,'%n%n%n%n');
            Um(1:length(freeS),j-1,i) = data{2}(freeS);
            Vm(1:length(freeS),j-1,i) = data{3}(freeS);
            Wm(1:length(freeS),j-1,i) = data{4}(freeS);
            if j == 3
                depthM(1:length(data{1}),i) = data{1};
            end
            fclose(fid);
            clear data
        end
    end
elseif avgLines == 1
    lineName = {'_alpha.water';'_U'};
    modelLegendText = {'Observations';'Averaged Lines (71.1 - 71.5)'};
    if length(RBRbursts) == 1
        burstTxt = sprintf('%d',RBRbursts(1));
    else
        burstTxt = sprintf('%d-%d',RBRbursts(1),RBRbursts(2));
    end
    modelBasePath = ['c:\Users\user\Documents\Models\MKK_Combined\2D\Burst' burstTxt '\' modelNames{1} '\MKKmodel\postProcessing\line\'];
    modelFolders = dir(modelBasePath);
    modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
    fileName = dir([modelBasePath modelFolders(1).name '\']);
    fileName = {fileName.name};
    isFile = contains(fileName,lineName{1});
    multiPlot = nnz(isFile);
    Em = NaN(profileLength,4096,multiPlot); %model epsilon
    Um = NaN(profileLength,4096,multiPlot); %model U
    Vm = NaN(profileLength,4096,multiPlot); %model V
    Wm = NaN(profileLength,4096,multiPlot); %model W
    depthM = NaN(profileLength,1,multiPlot);%model depth
    for i = 2:length(modelFolders)-1
        fileName = dir([modelBasePath modelFolders(i).name '\']);
        fileName = {fileName.name};
        %Load epsilon data
        isFile = find(contains(fileName,lineName{1}));
        for j = 1:multiPlot
            fid = fopen([modelBasePath modelFolders(i).name '\' fileName{isFile(j)}]);
            if contains(fileName{isFile(j)},'epsilon')
                data = textscan(fid,'%n%n%n%n%n');
                freeS = find(data{2} >= 0.5);
                %Determine if file is epsilon or omega (closure scheme model)
                Em(1:length(freeS),i-1,j) = data{3}(freeS);
            else
                data = textscan(fid,'%n%n%n%n%n%n');
                freeS = find(data{2} >= 0.5);
                %Determine if file is epsilon or omega (closure scheme model)
                k = data{3}(freeS);
                w = data{4}(freeS);
                Cmu = 0.09; %Constant
                Em(1:length(freeS),i-1,j) = Cmu.*k.*w;
            end
            
            fclose(fid);
            clear data
        end
        %Load U data
        isFile = find(contains(fileName,lineName{2}));
        for j = 1:multiPlot
            fid = fopen([modelBasePath modelFolders(i).name '\' fileName{isFile(j)}]);
            data = textscan(fid,'%n%n%n%n');
            Um(1:length(freeS),i-1,j) = data{2}(freeS);
            Vm(1:length(freeS),i-1,j) = data{3}(freeS);
            Wm(1:length(freeS),i-1,j) = data{4}(freeS);
            if i == 3
                depthM(1:length(data{1}),j) = data{1};
            end
            fclose(fid);
            clear data
        end
    end
    %Average lines together
    Um = nanmean(Um,3);
    Vm = nanmean(Vm,3);
    Wm = nanmean(Wm,3);
    multiPlot = 1;
elseif multipleModels == 1 && avgLines == 1
        error('Cannot plot multiple models and avg lines at the same time')
elseif multipleModels == 0 && avgLines == 0
        error('User must specify either mutlipleModels or avgLines')
end

%Compute wave statistics from model data
Uwrms_m = zeros(profileLength,multiPlot); %model wave velocities
Uw_m = zeros(profileLength,multiPlot);
Wwrms_m = zeros(profileLength,multiPlot);
Ww_m = zeros(profileLength,multiPlot);
E_m = zeros(profileLength,multiPlot);
for i = 1:multiPlot
    for j = 1:profileLength
        U = Um(j,:,i);
        V = Vm(j,:,i);
        W = Wm(j,:,i);
        
        %remove tailing NaNs
        U = U(~isnan(U));
        V = V(~isnan(V));
        W = W(~isnan(W));
        
        nwin = length(U);
        
        %Calculate RMS velocities (Lowe et al. 2005; Luhar et al. 2013)
        %Ec and Nc are mean east and north current velocities
        %Ewrms and Nwrms are rms oscillatory velocities
        %Uw and Ww are total mean oscillatory horizontal and vertical velocities
        Ec = (1/nwin)*sum(U);Nc = (1/nwin)*sum(V);Wc = (1/nwin)*sum(W);
        Ewrms = sqrt((1/nwin)*sum((U-Ec).^2));
        Nwrms = sqrt((1/nwin)*sum((V-Nc).^2));
        Wwrms_m(j,i) = sqrt((1/nwin)*sum((W-Wc).^2));
        Uwrms_m(j,i) = sqrt(Ewrms^2+Nwrms^2);
        Uw_m(j,i) = sqrt(2)*Uwrms_m(j,i);
        Ww_m(j,i) = sqrt(2)*Wwrms_m(j,i);
        E_m(j,i) = nanmean(Em(j,:,i));
    end
end
clear U V W
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

%Compute wave statistics per depth bin, windowed in time
avt = fs*512; %window 20 s long
idx = 1:avt:m*n;
Uwrms_i = zeros(length(data.depth),length(idx)-1); %observed wave velocities
Uw_i = zeros(length(data.depth),length(idx)-1);
Wwrms_i = zeros(length(data.depth),length(idx)-1);
Ww_i = zeros(length(data.depth),length(idx)-1);
E_i = zeros(length(data.depth),length(idx)-1);
for i = 1:length(idx)-1
    inds = idx(i):idx(i+1);
    for j = 1:length(data.depth)
        U = Ui(j,inds);
        V = Vi(j,inds);
        W = Wi(j,inds);
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
        Wwrms_i(j,i) = sqrt((1/nwin)*sum((W-Wc).^2));
        Uwrms_i(j,i) = sqrt(Ewrms^2+Nwrms^2);
        Uw_i(j,i) = sqrt(2)*Uwrms_i(j,i);
        Ww_i(j,i) = sqrt(2)*Wwrms_i(j,i);
    end
    
    %Compute dissipation rate of turbulence (Wiles et al., 2006)
    r = (cellSize)/cosd(25);   %along-beam distance (AQDPs have beam angles of 25 deg)
    lags = 5;
    [m,n] = size(Wi(:,inds));
    maxbin = m-lags;
    itt = 1;                   %iteration number
    D = zeros(maxbin,lags);    %structure function D(z,r)
    R = zeros(maxbin,lags);    %along-beam distance r^2/3
    epsilon = NaN(1,m); %dissipation rate [m^2/s^3]
    rr = zeros(1,lags);
    while itt <= maxbin
        d = zeros(n,lags);
        idl = itt:itt+lags;
        for kk = 1:lags %compute velocity differences (vel(z,1) - vel(z,2:lags))^2
            d(:,kk) = (Wi(itt,inds)-(Wi(idl(kk+1),inds))).^2;
            D(itt,kk) = nanmean(d(:,kk));
            rr(:,kk) = r*kk;
        end
        R(itt,:) = rr.^(2/3);
        [pf,S] = polyfit(R(itt,:),D(itt,:),1);  %linear regression of D along r^2/3
        A = pf(1);
        epsilon(idl(5)) = (A/2)^(3/2);
        itt = itt+1;
    end
    epsilon(abs(imag(epsilon))>0) = NaN;
    E_i(:,i) = epsilon;
end

%% Find depths consistent between the model and the instrument
cc = brewermap(multiPlot,'Set1');
p1 = zeros(1,1);
p2 = zeros(multiPlot,1);
% zuv(18:19) = NaN;
zuv = zuv-0.1;
depthM = depthM-0.08;

%Plot some simple comparisons for now
ff(1) = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[800 350   1250   500]);
for i = 1:multiPlot
    nearBed_m = 1:length(depthM(:,i));

    sp(1) = subplot(131);
    hold on
    if i == 1
        plot(nanmean(Uwrms_i,2),zuv,'+k','linewidth',1.5);
    end
    plot(Uwrms_m(nearBed_m,i),depthM(nearBed_m,i),'linewidth',1.5,'color',cc(i,:));
    sp(2) = subplot(132);
    hold on
    if i == 1
        plot(nanmean(Wwrms_i,2),zuv,'+k','linewidth',1.5);
    end
    plot(Wwrms_m(nearBed_m,i),depthM(nearBed_m,i),'linewidth',1.5,'color',cc(i,:));
    sp(3) = subplot(133);
    hold on
    if i == 1
        p1(1) = plot(nanmean(E_i,2),zuv,'+k','linewidth',1.5);
    end
    p2(i) = plot(E_m(nearBed_m,i),depthM(nearBed_m,i),'linewidth',1.5,'color',cc(i,:));
end

leg = legend([p1; p2],modelLegendText);
ylabel(sp(1),'\bf\itDepth [m]')
xlabel(sp(1),'\bf\itU_{w,RMS} [m/s]')
xlabel(sp(2),'\bf\itW_{w,RMS} [m/s]')
xlabel(sp(3),'\bf\it\epsilon [m^2/s^2]')
set(sp,'ylim',[0 0.7],'ytick',0:0.1:0.7)
set(sp(1),'xlim',[0 0.15])
set(sp(2),'xlim',[0 0.10])
set(sp(3),'xscale','log','xlim',[10^-7 10^-3])
set(sp(1),'position',[0.1 0.15 0.21 0.72])
set(sp(2),'position',[0.3805 0.15 0.21 0.72])
set(sp(3),'position',[0.65 0.15 0.21 0.72])
set(leg,'position',[0.88 0.6 0.05 0.05])
suptitle(['\bf\it' sprintf('Burst no. %s',burstTxt)])
prettyfigures('text',11,'labels',13)
if saveFigs == 1
    export_fig(ff(1),[figPath modelFigName '_AQDP_MODEL-Fields_V7'],'-png','-r600','-nocrop')
end
if saveData == 1
    dataPath = 'c:\Users\user\Documents\Models\DataAnalysis\Calibration\';
    clear data
    data.Gatts.modelName = modelNames;
    data.Gatts.instNames = insts;
    data.Gatts.lineName = lineName;
    data.Gatts.RBRbursts = RBRbursts;
    data.Gatts.AQDbursts = AQDbursts;
    data.depth_i = zuv;
    data.Uwrms_i = Uwrms_i;
    data.Wwrms_i = Wwrms_i;
    data.E_i = E_i;
    data.depth_m = depthM;
    data.Uwrms_m = Uwrms_m;
    data.Wwrms_m = Wwrms_m;
    data.E_m = E_m;
    save([dataPath modelFigName '_calibrationProfiles_V7.mat'],'data','-mat','-v7.3')
end