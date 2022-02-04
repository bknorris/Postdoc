%Compare pressure spectra from the model and selected instruments for
%validation of model functionality.
%
%
% This is version 7 of this script.
%
% Changelog:
% 12/21/20: using free surface derived from model to determine surface
% elevation spectrum instead of line sample
% 12/22/20: plots multiple instruments (instead of just AQDP and model)
% 01/06/21: plot multiple model results together to compare model settings
% 01/12/21: compute percent difference between observations and model
% results
% 01/29/21: updated script to combined multiple bursts togther to estimate
% spectra. Added switch to turn the pct diff text on the plots on or off.
% 03/02/21: updated free surface extraction method to extract only one z
% value per x value (takes mean of y axis). Double checked spectral calcs
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
lineName = 'alpha.water_freeSurface.raw';
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
modelNames = {'MKK_CombinedHR_V28'};
modelFigName = 'MKK_CombinedHR_V28';
insts = {'MKK18C1501rbr-b.nc';'MKK18HR101aqdHR-b.nc'};
instName = {'C15';'HR-Aquadopp'};
legendText = {'Observations';'V28, Level 2 @ free Surf'};

%Script Settings
InstxLoc = [29 71.3]; %C15, C16, HR AQDP
nfft1 = [128 256]; %spectral window settings for instruments
nfft2 = [320 320]; %spectral window settings for the model
RBRbursts = [151 158]; %RBR bursts to analyze (range: 141 - 162)
AQDbursts = [61 68]; %ADCP bursts to analyze (range: 51 - 72)
addPctError = 1; %add percent error text to plots; 0 is off, 1 is on
saveFigs = 0; %save figures to disk; 0 is off, 1 is on

%% 1. Extract free surface from model
freeS = NaN(length(modelNames),4096,length(InstxLoc));
for i = 1:length(modelNames)
    if length(RBRbursts) == 1
        burstTxt = sprintf('%d',RBRbursts(1));
    else
        burstTxt = sprintf('%d-%d',RBRbursts(1),RBRbursts(2));
    end
    modelBasePath = ['c:\Users\user\Documents\Models\MKK_Combined\2D\Burst' burstTxt '\' modelNames{i} '\MKKmodel\postProcessing\freeSurface\'];
    modelFolders = dir(modelBasePath);
    modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
    for j = 2:length(modelFolders)-1
        fileName = dir([modelBasePath modelFolders(j).name '\']);
        fileName = {fileName.name};
        isFile = find(contains(fileName,lineName));
        fid = fopen([modelBasePath modelFolders(j).name '\' fileName{isFile}]);
        data = textscan(fid,'%n%n%n%n','headerlines',2);
        [ud,ix,iy]=uniquetol(data{1},1e-4);
        output = [ud, accumarray(iy,data{3},[],@mean)];
        for k = 1:length(InstxLoc)
            [~,InstxID] = min(abs(output(:,1)-InstxLoc(k)));
            freeS(i,j,k) = output(InstxID,2);
        end
        fclose(fid);
        clear data
    end
end

%% 2. Compare model and instrument spectra
ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[500 350   950   500]);
sp = zeros(length(insts),1);
p1 = zeros(1,1);
p2 = zeros(length(modelNames),1);
cc = brewermap(length(modelNames),'Set1');
for i = 1:length(modelNames)
    FS = squeeze(freeS(i,1:512,:));
    FS = repmat(FS,8,1);
    FS = cmgbridge(FS,100,100,10000);

    %% Process instrument data
    for j = 1:length(insts)
        sp(j) = subplot(1,2,j);
        hold on;
        data = loadnc([dataPath insts{j}]);
        %Get attributes from the data
        if contains(data.Gatts.INST_TYPE,'RBR')
            lat = data.Gatts.latitude;                 %latitude for h calcs
            zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
            dt = 1/data.Gatts.sample_interval;         %sample interval [s]
            idx = find(data.burst >= RBRbursts(1) & data.burst <= RBRbursts(2));
            dof = (length(data.P_1)/nfft1(j))*2;
            fprintf('RBR spectra will be averaged with: %0.f Degrees of freedom\n',dof)
        else
            lat = data.Gatts.latitude;                 %latitude for h calcs
            zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
            dt = data.Gatts.instmeta_AQDSamplingrate;  %sample freq [Hz]
            idx = find(data.burst >= AQDbursts(1) & data.burst <= AQDbursts(2));
            dof = (length(data.P_1)/nfft1(j))*2;
            fprintf('ADCP spectra will be averaged with: %0.f Degrees of freedom\n',dof)
        end
        
        %Calculate mean water depth in burst
        [r,s] = size(data.dn(:,idx));
        P = reshape(data.P_1(:,idx),r*s,1)+zp;
        g = zeros(length(P),1);h = zeros(length(P),1);
        for l = 1:length(P)
            g(l,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*P(l,:);
            h(l,:) = ((((-1.82E-15*P(l,:)+2.279E-10)*P(l,:)-2.2512E-5)*P(l,:)+9.72659)*P(l,:))/g(l,:);
        end
        h = nanmean(h);                                 %mean depth during time record
        clear P r s
        fprintf('Mean water depth during time record: %0.2f m\n',h)
        
        %Spectral analysis
        [r,s] = size(data.dn(:,idx));
        P = zeros(r,s);
        for jj = 1:length(idx)
            P(:,jj) = detrend(data.P_1(:,idx(jj))); %detrend spectrum per segment
        end
        P = reshape(P,r*s,1);         %redefine for spectral calcs
        MP = detrend(FS(~isnan(FS(:,j)),j));
        [Cpp,F1] = pwelch(P,hanning(nfft1(j)),round(nfft1(j)*0.25,0),nfft1(j),dt);
        [Cid,F2] = pwelch(MP,hanning(nfft2(i)),round(nfft2(i)*0.25,0),nfft2(i),8);
        
        %Do surface elevation scaling for the instrument data
        df = F1(3)-F1(2);
        omega = 2*pi.*F1;
        k = qkhf(omega,h)./h;
        coshkhz = cosh(k*zp);
        coshkh = cosh(k*h);
        attn = coshkhz./coshkh;                     %transfer function for surface elevation scaling
        attn(attn<0.2) = 0.2;                       %limiter for transfer function
        Spp = Cpp./(attn.^2);                       %surface elevation spectrum
        lfc = find(F1>=0.033,1,'first');hfc = find(F1<=0.33,1,'last'); %lo and hi freq cutoffs for Hs calcs
        HsI = 4*sqrt(sum(Spp(lfc:hfc).*mean(diff(F1(lfc:hfc))))); %Hs from the instrument
        lfc = find(F2>=0.033,1,'first');hfc = find(F2<=0.33,1,'last');
        HsM = 4*sqrt(sum(Cid(lfc:hfc).*mean(diff(F2(lfc:hfc))))); %Hs from the model
        pctErr = ((HsM-HsI)/HsI)*100;
                
        if i == 1
            p1(i) = plot(F1,Spp,'+k','linewidth',1.5);
            if addPctError == 1
                text(0.2,0.035,sprintf('Pct. Error: %0.2f',pctErr))
            end
        end
        p2(i) = plot(F2,Cid,'color',cc(i,:),'linewidth',1.5);

    end
end
leg = legend([p1; p2],legendText);
set(leg,'position',[0.7 0.7 0.25 0.1])
ylabel(sp(1),'\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
xlabel(sp(1),'\bf\itf [Hz]'),xlabel(sp(2),'\bf\itf [Hz]')
set(sp,'xlim',[0 0.33],'ylim',[0 0.06])
title(sp(1),['\bf\itModel Inlet - ' instName{1}])
title(sp(2),['\bf\itHi-Res Patch - ' instName{2}])

ff2 = figure(2);
set(ff2,'PaperOrientation','landscape',...
    'position',[100 80   650   400]);
p = zeros(length(modelNames),1);
hold on
for i = 1:length(modelNames)
    xx = linspace(1,(length(freeS)/8),length(freeS));
    p(i) = plot(xx,freeS(i,:,2),'color',cc(i,:),'linewidth',1.5);
end
ylabel('\bf\itWater Depth [m]')
xlabel('\bf\itTime [s]')
set(gca,'ylim',[-0.35 0.25],'xlim',[0 64],'position',[0.1 0.15 0.8 0.75])
title('\bf\itModeled Free Surface')
leg = legend(p,legendText{2:end});
set(leg,'position',[0.75 0.25 0.1 0.1])

prettyfigures('text',11,'labels',12,'box',1)
if saveFigs == 1
    export_fig(ff1,[figPath  modelFigName 'burst' burstTxt '_AQDP_MODEL-spectra_V7'],'-png','-r600','-nocrop')
    export_fig(ff2,[figPath  modelFigName 'burst' burstTxt '_AQDP_MODEL-freeSurf_V7'],'-png','-r600','-nocrop')
end









