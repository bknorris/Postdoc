%Compare pressure spectra from the model and selected instruments for
%validation of model functionality.
%
%
% This is version 5 of this script.
%
% Changelog:
% 12/21/20: using free surface derived from model to determine surface
% elevation spectrum instead of line sample
% 12/22/20: plots multiple instruments (instead of just AQDP and model)
% 01/06/21: plot multiple model results together to compare model settings
% 01/12/21: compute percent difference between observations and model
% results
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
lineName = 'alpha.water_freeSurface.raw';
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
modelFigName = 'MKK_CombinedHR_V10';
modelNames = {'MKK_CombinedHR_V10'};
insts = {'MKK18C1501rbr-b.nc';'MKK18HR101aqdHR-b.nc'};
InstxLoc = [29 69.5]; %C15, C16, HR AQDP
whichBurst = [160 70]; %burst to perform analysis on (hint: check MKK_ModelLogs.xlsx)

%% Extract data from the model directory
freeS = zeros(length(modelNames),4096,length(InstxLoc));
for i = 1:length(modelNames)
    modelBasePath = ['c:\Users\user\Documents\Models\MKK_Combined\2D\Burst' num2str(whichBurst(1)) '\' modelNames{i} '\MKKmodel\postProcessing\freeSurface\'];
    modelFolders = dir(modelBasePath);
    modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
    for j = 2:length(modelFolders)-1
        fileName = dir([modelBasePath modelFolders(j).name '\']);
        fileName = {fileName.name};
        isFile = find(contains(fileName,lineName));
        fid = fopen([modelBasePath modelFolders(j).name '\' fileName{isFile}]);
        data = textscan(fid,'%n%n%n%n','headerlines',2);
        for k = 1:length(InstxLoc)
            [~,InstxID] = min(abs(data{1}-InstxLoc(k)));
            freeS(i,j,k) = data{3}(InstxID);
        end
        fclose(fid);
        clear data
    end
end

%% Set up the figure
ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[500 350   950   500]);
sp = zeros(length(insts),1);
p = zeros(length(modelNames),1);
m = zeros(1,1);
cc = brewermap(4,'Set1');
for i = 1:length(modelNames)
    FS = squeeze(freeS(i,:,:));
    %% Now load the instrument data
    for j = 1:length(insts)
        sp(j) = subplot(1,2,j);
        hold on;
        data = loadnc([dataPath insts{j}]);
        burstID = find(data.burst == whichBurst(j));
        %Load the data, get attributes
        if contains(data.Gatts.INST_TYPE,'RBR')
            lat = data.Gatts.latitude;                 %latitude for h calcs
            zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
            dt = 1/data.Gatts.sample_interval;         %sample rate [Hz]
        else
            lat = data.Gatts.latitude;                 %latitude for h calcs
            zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
            dt = data.Gatts.instmeta_AQDSamplingrate;  %sample rate [Hz]
        end
        
        %Calculate mean water depth in burst
        time = data.dn(:,burstID);
        P = data.P_1(:,burstID);
        g = zeros(length(P),1);h = zeros(length(P),1);
        for l = 1:length(P)
            g(l,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*P(l,:);
            h(l,:) = ((((-1.82E-15*P(l,:)+2.279E-10)*P(l,:)-2.2512E-5)*P(l,:)+9.72659)*P(l,:))/g(l,:);
        end
        h = mean(h);
        
        %% Compute spectra
        P = detrend(P);
        MP = detrend(FS(:,j));
        if j == 1
            nfft1 = 128;
            nfft2 = 512;
        else
            nfft1 = 128;
            nfft2 = 300;
        end
        [Cpp,F1] = pwelch(P,hanning(nfft1),round(nfft1*0.7,0),nfft1,dt);
        [Cid,F2] = pwelch(MP,hanning(nfft2),round(nfft2*0.7,0),nfft2,8);
        
        %Do surface elevation scaling for the instrument data
        df = F1(3)-F1(2);
        omega = 2*pi.*F1;
        k = qkhf(omega,h)./h;
        coshkhz = cosh(k*zp);
        coshkh = cosh(k*h);
        attn = coshkhz./coshkh;                     %transfer function for surface elevation scaling
        attn(attn<0.2) = 0.2;                       %limiter for transfer function
        Spp = Cpp./(attn.^2);                       %surface elevation spectrum
        HsI = 4*sqrt(sum(Spp.*mean(diff(F1)))); %Hs from the instrument
        HsM = 4*sqrt(sum(Cid.*mean(diff(F2)))); %Hs from the model
        pctErr = ((HsM-HsI)/HsI)*100;
        
        if i == 1
            m(i) = plot(F1,Spp,'+k','linewidth',1.5);
            text(0.2,0.015,sprintf('Pct. Error: %0.2f',pctErr))
        end
        p(i) = plot(F2,Cid,'color',cc(i,:),'linewidth',1.5);
    end
end
leg = legend([m; p],{'Observations';'Model'});
set(leg,'position',[0.7 0.7 0.25 0.1])
ylabel(sp(1),'\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
xlabel(sp(1),'\bf\itf [Hz]'),xlabel(sp(2),'\bf\itf [Hz]')
set(sp,'xlim',[0 0.5],'ylim',[0 0.1])
title(sp(1),insts{1})
title(sp(2),insts{2})

prettyfigures('text',11,'labels',12,'box',1)
% export_fig(ff1,[figPath  modelFigName 'burst' num2str(whichBurst(1)) '_AQDP_MODEL-spectra_V5'],'-png','-r600','-nocrop')

