%Compare pressure spectra from the model and selected instruments for
%validation of model functionality.
%
%
% This is version 3 of this script.
%
% Changelog:
% 12/21/20: using free surface derived from model to determine surface
% elevation spectrum instead of line sample
% 12/22/20: plots multiple instruments (instead of just AQDP and model)
%
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
modelBasePath = 'c:\Users\user\Documents\Models\MKK_Combined\2D\MKK_CombinedHR_V8\MKKmodel\postProcessing\freeSurface\';
% modelBasePath = 'c:\Users\user\Documents\Models\MKK_Combined\2D\MKK_CombinedLR_V1\MKKmodel\postProcessing\freeSurface\';
lineName = 'alpha.water_freeSurface.raw';
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
insts = {'MKK18C1501rbr-b.nc';'MKK18HR101aqdHR-b.nc'};
% insts = {'MKK18C1501rbr-b.nc';'MKK18LR101aqdHR-b.nc'};
InstxLoc = [19 71.5]; %C15, C16, HR AQDP
% InstxLoc = [29 69]; %C15, C16, LR AQDP
whichBurst = [144 54]; %burst to perform analysis on (hint: check MKK_ModelLogs.xlsx)
modelName = 'MKK_CombinedHR_V8';
% modelName = 'MKK_CombinedLR_V1';

%% Extract data from the model directory
modelFolders = dir(modelBasePath);
modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
freeS = zeros(length(modelFolders)-1,length(InstxLoc));
for i = 2:length(modelFolders)-1
    fileName = dir([modelBasePath modelFolders(i).name '\']);
    fileName = {fileName.name};
    isFile = find(contains(fileName,lineName));
    fid = fopen([modelBasePath modelFolders(i).name '\' fileName{isFile}]);
    data = textscan(fid,'%n%n%n%n','headerlines',2);
    for j = 1:length(InstxLoc)
        [~,InstxID] = min(abs(data{1}-InstxLoc(j)));
        freeS(i,j) = data{3}(InstxID);
    end
    fclose(fid);
    clear data
end

%% Set up the figure
ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[500 350   950   500]);
sp = zeros(length(insts),1);
%% Now load the instrument data
for i = 1:length(insts)
    data = loadnc([dataPath insts{i}]);
    burstID = find(data.burst == whichBurst(i));
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
    MP = detrend(freeS(:,i));
    if i == 1
        nfft1 = 128;
        nfft2 = 512;
    else
        nfft1 = 256;
        nfft2 = 256;
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
    hsI = 4*sqrt(sum(Spp.*mean(diff(F1)))); %Hs from the instrument
    hsM = 4*sqrt(sum(Cid.*mean(diff(F2)))); %Hs from the model
    
    sp(i) = subplot(1,2,i);
    p(1) = plot(F1,Spp,'b','linewidth',1.5);
    hold on
    p(2) = plot(F2,Cid,'r','linewidth',1.5);
    text(0.2,0.05,sprintf('H_s (observations): %0.2f',hsI))
    text(0.2,0.043,sprintf('H_s (model): %0.2f',hsM))
    leg = legend(p,{'Instrument';'Model'});
    set(leg,'box','off')
    hold off;
end
ylabel(sp(1),'\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
xlabel(sp(1),'\bf\itf [Hz]'),xlabel(sp(2),'\bf\itf [Hz]')
set(sp,'xlim',[0 0.5],'ylim',[0 0.1])
title(sp(1),insts{1})
title(sp(2),insts{2})

ff2 = figure(2);
    set(ff2,'PaperOrientation','landscape',...
    'position',[100 80   650   400]);
xx = linspace(1,(length(freeS)/8),length(freeS));
plot(xx,freeS(:,2),'b','linewidth',1.5)
ylabel('\bf\itWater Depth [m]')
xlabel('\bf\itTime [s]')
set(gca,'ylim',[-0.25 0.25],'xlim',[3 512],'position',[0.1 0.15 0.8 0.75])
suptitle('Modeled Free Surface')

prettyfigures('text',11,'labels',12,'box',1)
export_fig(ff1,[figPath  modelName '_AQDP_MODEL-spectra_V3'],'-png')
export_fig(ff2,[figPath modelName '_MODEL-freeSurf'],'-png')

