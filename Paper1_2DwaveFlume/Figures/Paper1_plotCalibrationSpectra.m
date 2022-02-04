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
modelPath = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\Calibration\';
lineName = 'alpha.water_freeSurface.raw';
figPath = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\Figures\'; %where to save figs
modelNames = {'MKK_CombinedHR_V27';'MKK_CombinedLR_V2'};
modelFigName = 'MKK_CombinedHR_V27';
insts = {'MKK18HR101aqdHR-b.nc';'MKK18LR101aqdHR-b.nc'};
instName = {'HR \ Site';'LR \ Site'};
legendText = {'Observation';'Model'};

%Script Settings
InstxLoc = [71.3 67]; %HR AQDP, LR AQDP
nfft1 = [128 128]; %spectral window settings for instruments
nfft2 = [320 320]; %spectral window settings for the model
RBRbursts = [151 158]; %RBR bursts to analyze (range: 141 - 162)
AQDbursts = [61 68]; %ADCP bursts to analyze (range: 51 - 72)
addPctError = 1; %add percent error text to plots; 0 is off, 1 is on
saveFigs = 0; %save figures to disk; 0 is off, 1 is on

%% 1. Extract free surface from model
freeS = NaN(length(modelNames),4096);
for i = 1:length(modelNames)
    modelBasePath = [modelPath '\' modelNames{i} '\MKKmodel\postProcessing\freeSurface\'];
    modelFolders = dir(modelBasePath);
    modelFolders = modelFolders(~ismember({modelFolders.name},{'.','..'}));
    modelFolders = {modelFolders.name};
    
    %Sort modelFolders
    sortModels = str2double(modelFolders); %windows sorts 1 then 10 then 100 then 2 then 20, etc.
    [~,sortID] = sort(sortModels);
    
    for j = 1:length(sortID)
        fileName = dir([modelBasePath modelFolders{sortID(j)} '\']);
        fileName = {fileName.name};
        isFile = find(contains(fileName,lineName));
        fid = fopen([modelBasePath modelFolders{sortID(j)} '\' fileName{isFile}]);
        data = textscan(fid,'%n%n%n%n','headerlines',2);
        [ud,ix,iy]=uniquetol(data{1},1e-4);
        output = [ud, accumarray(iy,data{3},[],@mean)];
        [~,InstxID] = min(abs(output(:,1)-InstxLoc(i)));
        freeS(i,j) = output(InstxID,2);
        fclose(fid);
        clear data
    end
end

%% 2. Compare model and instrument spectra
ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[500 350   650   350]);
sp = zeros(length(insts),1);
p1 = zeros(1,1);
p2 = zeros(length(modelNames),1);
cc = brewermap(1,'Set1');
marker = {'+';'x'};
for i = 1:length(modelNames)
    FS = freeS(i,1:512);
    FS = repmat(FS,1,8);
    FS = cmgbridge(FS,100,100,10000);
    
    %% Process instrument data
    sp(i) = subplot(1,2,i);
    hold on;
    data = loadnc([dataPath insts{i}]);
    %Get attributes from the data
    lat = data.Gatts.latitude;                 %latitude for h calcs
    zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
    dt = data.Gatts.instmeta_AQDSamplingrate;  %sample freq [Hz]
    idx = find(data.burst >= AQDbursts(1) & data.burst <= AQDbursts(2));
    dof = (length(data.P_1)/nfft1(i))*2;
    fprintf('ADCP spectra will be averaged with: %0.f Degrees of freedom\n',dof)
    
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
    MP = detrend(FS(~isnan(FS)));
    [Cpp,F1] = pwelch(P,hanning(nfft1(i)),round(nfft1(i)*0.25,0),nfft1(i),dt);
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
    
    p1(i) = plot(F1,Spp,marker{i},'markeredgecolor','k','linewidth',1.5);
    text(0.18,0.038,['$\mathrm{' sprintf('%s',instName{i}) '}$'],'fontsize',10,'interpreter','latex')
    text(0.18,0.035,['$\mathrm{' sprintf('Error: %0.2f',pctErr) '}\%$'],'fontsize',8,'interpreter','latex')
    p2(i) = plot(F2,smooth(Cid,2),'color',cc(1,:),'linewidth',1.5);

end
%Legend
leg = legend([p1(1); p2(1)],legendText);
set(leg,'position',[0.25 0.8 0.25 0.1],'box','off')

%Labeling
ylabel(sp(1),'$S(f)_{\eta\eta} \quad \mathrm{(m^2/Hz)}$','interpreter','latex')
xlabel(sp(1),'$f\ \mathrm{(Hz)}$','interpreter','latex')
xlabel(sp(2),'$f\ \mathrm{(Hz)}$','interpreter','latex')

%Positioning
set(sp,'xlim',[0 0.33],'ylim',[0 0.06])
set(sp(2),'yticklabel',[])
set(sp(1),'position',[0.1 0.15 0.38 0.78])
set(sp(2),'position',[0.52 0.15 0.38 0.78])

prettyfigures('text',11,'labels',10)
set(ff1,'units','inches');
pos = get(ff1,'Position');
set(ff1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% export_fig(ff1,[figPath modelFigName '_AQDP_MODEL-Cals_V2'],'-pdf','-nocrop','-nofontswap')
print(ff1,[figPath modelFigName '_AQDP_MODEL-SpectraCals_V2'],'-dpdf','-r0')









