%Compare pressure spectra from the model and selected instruments for
%validation of model functionality.
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
modelBasePath = 'c:\Users\user\Documents\Models\MKK_CombinedHR\2D\MKK_CombinedLR_V1\MKKmodel\postProcessing\line\';
lineName = 'line2_alpha.water_epsilon_p_rgh.xy';
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
insts = 'MKK18LR101aqdHR-b.nc';
modelName = 'MKK_CombinedLR_V1';

%% Extract data from the model directory
rampTime = 128; %find values greater than rampTime
modelFolders = dir(modelBasePath);
idx = cellfun(@(x)str2double(x) > rampTime,{modelFolders.name},'UniformOutput',false);
idx = cell2mat(idx);
idx2 = find(idx);
freeS = zeros(1,length(modelFolders)-3);
for i = 4:length(modelFolders) %skip . and ..
    fileName = dir([modelBasePath modelFolders(i).name '\']);
    fileName = {fileName.name};
    isFile = find(contains(fileName,lineName));
    fid = fopen([modelBasePath modelFolders(i).name '\' fileName{isFile}]);
    data = textscan(fid,'%n%n%n%n');
    surfID = find(data{2} > 0.001,1,'last');
    freeS(i) = data{1}(surfID);
    fclose(fid);
    clear data
end

%% Now load the instrument data
data = loadnc([dataPath insts]);
whichBurst = 54; %burst to perform analysis on (hint: check MKK_ModelLogs.xlsx)
burstID = find(data.burst == whichBurst);
%Load the data, get attributes

lat = data.Gatts.latitude;                 %latitude for h calcs
zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
dt = data.Gatts.instmeta_AQDSamplingrate;  %sample rate [Hz]

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
% MP = detrend(freeS(idx2));
MP = detrend(freeS(idx2));
[Cpp,F1] = pwelch(P,hanning(round((length(P)/36))),round((length(P)/36)*0.7,0),round((length(P)/36),0),dt);
[Cid,F2] = pwelch(MP,hanning(round((length(MP)/12))),round((length(MP)/12)*0.7,0),round((length(MP)/12),0),8);

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

ff(1) = figure(1);
    set(ff,'PaperOrientation','landscape',...
    'position',[500 350   550   500]);
p(1) = plot(F1,Spp,'b','linewidth',1.5);
hold on
p(2) = plot(F2,Cid,'r','linewidth',1.5);
text(0.35,0.05,sprintf('H_s (ADCP): %0.2f',hsI))
text(0.35,0.043,sprintf('H_s (model): %0.2f',hsM))
ylabel('\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
xlabel('\bf\itf [Hz]')
set(gca,'xlim',[0 0.5])
leg = legend(p,{'HR Aquadopp';'Model'});
set(leg,'box','off')
hold off;


ff(2) = figure(2);
    set(ff,'PaperOrientation','landscape',...
    'position',[100 80   650   400]);
xx = linspace(1,(length(freeS)/8),length(freeS));
p(1) = plot(xx,freeS,'b','linewidth',1.5);
ylabel('\bf\itWater Depth [m]')
xlabel('\bf\itTime [s]')
set(gca,'ylim',[2 3],'xlim',[3 512],'position',[0.1 0.15 0.8 0.75])
suptitle('Modeled Free Surface')
prettyfigures('text',11,'labels',12,'box',1)

export_fig(ff(1),[figPath  modelName '_AQDP_MODEL-spectra'],'-png')
export_fig(ff(2),[figPath modelName '_MODEL-freeSurf'],'-png')

%
