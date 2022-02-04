%% Calculate Wave Inputs for IHFOAM irregularWaveModel
% The model requires wave period, height (amplitude*2), phase (radians),
% and direction
%
% NOTE: C15 doesn't measure (u,v), so we will assume the wave field between
%       the instrument and the hi-res model is shore-perpendicular.
%
% Inputs: pressure time series from C15
% Outputs: formatted waveProperties for the model and optional plot
%
% Updates: 12/01/2020: fixed definition for Hi
%          01/06/2021: added a limter for frequencies < 0.5 Hz to reduce
%          number of output frequencies. Optimal nfreq ~ 200 
% This is version 2 of this script.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% close all;

%% 1. Plot wave statistics from C14 to determine which burst to analyze
dataPath = 'C:\Users\user\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs

load([dataPath 'MKK18C1501rbr_waves.mat'])
ploton = 0; %turns plots on (= 1) or off (= 0)

%% 2. Pick a Burst to Analyze
burstNum = 160;

%Load the data, get attributes
data = loadnc([dataPath 'MKK18C1501rbr-b.nc']);
lat = data.Gatts.latitude;                 %latitude for h calcs
zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
dt = data.Gatts.sample_interval;           %sample interval [s]
burstID = find(data.burst == burstNum);

%Calculate mean water depth in burst
time = data.dn(:,burstID);
P = data.P_1(:,burstID)+zp;
g = zeros(length(P),1);h = zeros(length(P),1);
for l = 1:length(P)
    g(l,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*P(l,:);
    h(l,:) = ((((-1.82E-15*P(l,:)+2.279E-10)*P(l,:)-2.2512E-5)*P(l,:)+9.72659)*P(l,:))/g(l,:);
end
h = nanmean(h);                     %mean depth during burst

%% 3. Spectral analysis to derive T_i, H_i, and psi_i

%Settings
m = 640;                          %number of windows to average spectra, max(T) output will be m/2
noverlap = 320;                   %number of overlap points
P = data.P_1(:,burstID);          %redefine P for spectral calcs
%FFT
n = length(P);
Fs = 1./dt;
dof = (n/m)*2;
fprintf('Spectra will be averaged with: %0.f Degrees of freedom\n',dof)

%%% The following lines from calcspectrum.m %%%
index = 1:m;
kk = fix((n-noverlap)/(m-noverlap));	% Number of windows
w = hanning(m);                     %Window specification
KMU = kk*norm(w)^2;                  %Normalizing scale factor
Pxx = zeros(m,1);
Phx = zeros(m,1);
%Compute windowed spectra
for i=1:kk
    xw = w.*detrend(P(index));
    index = index + (m - noverlap);
    Xx = abs(fft(xw)).^2;
    Ph = angle(fft(xw));
    Pxx = Pxx + Xx;
    Phx = Phx + Ph;
end
% Select first half and let 1st point (DC value) be replaced
% with 2nd point
select = [2 2:m/2];
Pxx = Pxx(select);
Phx = Phx(select);
Cpp = Pxx./KMU;
f = (select - 1)'*Fs/m;

%%% The following lines from waveStats or wds.m %%%
%Depth correction for pressure sensor elevation
omega = 2*pi.*f;
k = qkhf(omega,h)./h;                %wavenumber from depth
coshkhz = cosh(k*zp);
coshkh = cosh(k*h);
attn = coshkhz./coshkh;              %transfer function for surface elevation scaling
attn(attn<0.2) = 0.2;                %limiter for transfer function
Spp = Cpp./(attn.^2);                %surface elevation spectrum
fid = find(f >= 0.0156 & f <= 0.5);  %limit output to f < 0.5 Hz 
f = f(fid);
Spp = Spp(fid);
Phx = Phx(fid); 

%%% The following lines adapted from ubspec.m %%%
Hs = 4*sqrt(sum(Spp.*mean(diff(f))));  %To compare with waveStats
pctDiff = ((wave.Hs(burstID)-Hs)/(Hs))*100;
fprintf('Percent difference in Hs between waveStats and this script: %0.2f %%\n\n',pctDiff)

%%% Quantities required for model inputs %%%
amp = sqrt(2*Spp.*mean(diff(f)));     %Wave amp = sqrt(2*S(fi)*df)
Hi = 2.*amp;                          %RMS Wave height = 2*amp 
Ti = 1./f;                            %Wave period per freq
Psi = wrapTo2Pi(Phx);                 %Wave phase between [0 2pi]

%Print out formatted data blocks to copy/paste into waveProperties
fprintf('waveHeights\n%0.0f(\n',length(Hi))
fprintf('(%0.6f)\n',Hi)
fprintf(');\n'), pause, clc

fprintf('wavePeriods\n%0.0f(\n',length(Ti))
fprintf('(%0.6f)\n',Ti)
fprintf(');\n'), pause, clc

fprintf('wavePhases\n%0.0f(\n',length(Psi))
fprintf('(%0.6f)\n',Psi)
fprintf(');\n'), pause, clc

fprintf('waveDirs\n%0.0f(\n',length(Hi))
fprintf('(%0.0f)\n',zeros(length(Hi),1))
fprintf(');\n'), pause, clc, disp('finished!')

%Plot Routine
if ploton == 1
    f1 = figure(1);
    set(f1,'PaperOrientation','landscape',...
        'position',[100 80   850   550]);
    tidx = find(wave.time>=wave.atts.AQDPstart&wave.time<=wave.atts.AQDPend);
    cc = brewermap(3,'Paired');
    sp(1) = subplot(3,1,1);
    plot(wave.burst(tidx),wave.depth(tidx),'-k','linewidth',1.5)
    sp(2) = subplot(3,1,2);
    p(1) = plot(wave.time(tidx),wave.Hs(tidx),'color',cc(1,:),'linewidth',1.5);hold on
    p(2) = plot(wave.time(tidx),wave.HrmsSS(tidx),'color',cc(2,:),'linewidth',1.5);
    p(3) = plot(wave.time(tidx),wave.HrmsIG(tidx),'color',cc(3,:),'linewidth',1.5);
    leg = legend(p,{'H_s';'H_{rms,SS}';'H_{rms,IG}'});
    sp(3) = subplot(3,1,3);
    p(1) = plot(wave.time(tidx),wave.Tp(tidx),'color',cc(1,:),'linewidth',1.5);hold on
    p(2) = plot(wave.time(tidx),wave.Tm(tidx),'color',cc(2,:),'linewidth',1.5);
    p(3) = plot(wave.time(tidx),wave.Tz(tidx),'color',cc(3,:),'linewidth',1.5);
    
    set(leg,'position',[0.88 0.49 0.05 0.05])
    ylabel(sp(1),'\bf\it Depth [m]')
    ylabel(sp(2),'\bf\it Wave Height [m]')
    ylabel(sp(3),'\bf\it Wave Period [s]')
    xlabel(sp(1),'\bf\itBurst Number')
    datetickzoom(sp(3),'x','dd-mm HH:MM:SS')
    set(sp(1),'ylim',[0 2],'ytick',0:0.5:2,'position',[0.1 0.75 0.75 0.2],'xlim',[min(wave.burst(tidx)) max(wave.burst(tidx))])
    set(sp(2),'ylim',[0 0.6],'ytick',0:0.2:0.6,'position',[0.1 0.4 0.75 0.2],'xlim',[min(wave.time(tidx)) max(wave.time(tidx))],'xticklabel',[])
    set(sp(3),'position',[0.1 0.12 0.75 0.2],'xlim',[min(wave.time(tidx)) max(wave.time(tidx))])
    set(leg,'position',[0.9 0.33 0.05 0.05])
    
    f2 = figure(2);
    set(f2,'PaperOrientation','landscape',...
        'position',[100 80   500   450]);
    p(1) = plot(f,Spp,'b','linewidth',1.5);
    ylabel('\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
    xlabel('\bf\itf [Hz]')
    set(gca,'xlim',[0 0.5])
    title(['\bf\itC15 Spectra - burst no. ' num2str(burstNum)])
    prettyfigures('text',12,'labels',13,'box',1,'tickdir','out','tlength',[0.005 0.005])
    export_fig(f1,[figPath 'MKK_C15_burst' num2str(burstNum) '_Conditions'],'-png')
    export_fig(f2,[figPath 'MKK_C15_burst' num2str(burstNum) '_spectra'],'-png')
end
