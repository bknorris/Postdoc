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
%          01/27/2021: redesign of script to plot figs first, compute input
%          spectra from multiple bursts for a more 'representative'
%          selection of time
%          03/02/2021: fixed spectral calcs to match output of pwelch
%
% This is version 4 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all;

%% 1. Plot wave statistics from C15 and the Aquadopp to determine which burst to analyze
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\user\Documents\Data\InstrumentData\Figures\'; %where to save figs
insts = {'MKK18C1501rbr-b.nc';'MKK18HR101aqdHR-b.nc'}; %'MKK18HR101aqdHR-b.nc'
instName = {'C15';'HR-Aquadopp'};
startTime = datenum('24-Jun-2018 21:00:00'); %RBR start time
endTime = datenum('25-Jun-2018 18:34:08'); %ADCP stop time
nfft = [128 256]; %spectral window settings
RBRbursts = [151 158]; %RBR bursts to analyze (range: 141 - 162)
AQDbursts = [61 68]; %ADCP bursts to analyze (range: 51 - 72)
printModelInput = 0; %display model input; 0 turns off; 1 turns on
saveFigs = 0; %save figures; 0 turns off; 1 turns on

%Create overview figure
ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[600 350   550   600],...
    'color','w','paperpositionmode','auto');
sp = zeros(3,1);pp = zeros(2,1);
cc1 = [29 145 192;34 94 168]./255;
cc2 = brewermap(300,'*YlGnBu');
for i = 1:length(insts)
    data = loadnc([dataPath insts{i}]);
    %Get attributes from the data
    if contains(data.Gatts.INST_TYPE,'RBR')
        lat = data.Gatts.latitude;                 %latitude for h calcs
        zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
        Fs = 1/data.Gatts.sample_interval;         %sample interval [s]
        idx = find(data.burst >= RBRbursts(1) & data.burst <= RBRbursts(2));
    else
        lat = data.Gatts.latitude;                 %latitude for h calcs
        zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
        Fs = data.Gatts.instmeta_AQDSamplingrate;  %sample rate [Hz]
        idx = find(data.burst >= AQDbursts(1) & data.burst <= AQDbursts(2));
    end
    %Calculate mean water depth in burst
    time = zeros(length(data.burst),1);
    depth = zeros(length(data.burst),1);
    Spp = zeros(nfft(i)/2+1,length(data.burst));
    lf = 0.03;hf = 0.3; %low and high freq cutoff values [Hz]
    
    for j = 1:length(data.burst)
        %Calculate mean water depth in burst
        time(j) = median(data.dn(:,j));
        P = data.P_1(:,j);
        P = cmgbridge(P,100,1000,10000);
        g = zeros(length(P),1);h = zeros(length(P),1);
        for l = 1:length(P)
            g(l,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*P(l,:);
            h(l,:) = ((((-1.82E-15*P(l,:)+2.279E-10)*P(l,:)-2.2512E-5)*P(l,:)+9.72659)*P(l,:))/g(l,:);
        end
        depth(j) = mean(h);
        h = mean(h);
        
        %Compute spectra using pwelch for efficiency
        P = detrend(P);
        [Cpp,F] = pwelch(P,hanning(nfft(i)),round(nfft(i)*0.25,0),nfft(i),Fs);
        
        %Surface elevation scaling for instrument data
        lfc = find(F >= lf,1,'first');              %low freq cutoff
        hfc = find(F <= hf,1,'last');               %high freq cutoff
        df = F(3)-F(2);
        omega = 2*pi.*F;
        k = qkhf(omega,h)./h;
        coshkhz = cosh(k*zp);
        coshkh = cosh(k*h);
        attn = coshkhz./coshkh;                     %transfer function for surface elevation scaling
        attn(attn<0.2) = 0.2;                       %limiter for transfer function
        Spp(:,j) = Cpp./(attn.^2);                  %surface elevation spectrum
    end
    
    tidx = find(time >= startTime & time <= endTime);
    sp(1) = subplot(3,1,1);
    pp(i) = plot(time(tidx),depth(tidx),'linewidth',1.5,'color',cc1(i,:));
    hold on
    sp(i+1) = subplot(3,1,i+1);
    imagesc(data.burst(tidx),F,Spp(:,tidx))
    hold on
    patch([data.burst(idx(1)) data.burst(idx(end)) data.burst(idx(end)) data.burst(idx(1)) data.burst(idx(1))],...
        [0.33 0.33 0 0 0.33],'r','facealpha',0.15)
    plot([data.burst(idx(1)) data.burst(idx(end)) data.burst(idx(end)) data.burst(idx(1)) data.burst(idx(1))],...
        [0.33 0.33 0 0 0.33],'color','r','linewidth',1.5)
    set(gca,'ylim',[0 0.33],...
        'xlim',[min(data.burst(tidx)) max(data.burst(tidx))])
    caxis([0.01 0.1])
    colormap(cc2);
    title(['\bf\it' instName{i}])
    clear data
end
%Legend
leg = legend(pp,instName,'location','eastoutside');
set(leg,'position',[0.82 0.91 0.05 0.02])

%Labeling
datetick(sp(1),'x','mm-dd HH:MM')
xlabel(sp(1),'\bf\itDateTime in 2018')
xlabel(sp(2),'\bf\itBurst no.')
xlabel(sp(3),'\bf\itBurst no.')
ylabel(sp(1),'\bf\itDepth (m)')
ylabel(sp(2),'\bf\itf (Hz)')
ylabel(sp(3),'\bf\itf (Hz)')

%Colorbar
cb = colorbar;
ylabel(cb,'\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
set(cb,'position',[0.82 0.13 0.025 0.5],'linewidth',1.5,'tickdir','out')

%Positioning, adjustments
set(sp(1),'position',[0.1 0.82 0.7 0.14],...
    'xlim',[startTime endTime],...
    'ylim',[0.75 2],'ytick',0.75:0.5:2)
set(sp(2),'position',[0.1 0.45 0.7 0.25])
set(sp(3),'position',[0.1 0.08 0.7 0.25])

%% 2. Compute model inputs from spectra
ff2 = figure(2);
set(ff2,'PaperOrientation','landscape',...
    'position',[1200 350   450   400],...
    'color','w','paperpositionmode','auto');
for i = 1:length(insts)
    data = loadnc([dataPath insts{i}]); %Load the data
    %Get attributes from the data
    if contains(data.Gatts.INST_TYPE,'RBR')
        lat = data.Gatts.latitude;                 %latitude for h calcs
        zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
        Fs = 1/data.Gatts.sample_interval;         %sample interval [s]
        idx = find(data.burst >= RBRbursts(1) & data.burst <= RBRbursts(2));
        dof = (length(data.P_1)/nfft(i))*2;
        fprintf('RBR spectra will be averaged with: %0.f Degrees of freedom\n',dof)
    else
        lat = data.Gatts.latitude;                 %latitude for h calcs
        zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
        Fs = data.Gatts.instmeta_AQDSamplingrate;  %sample freq [Hz]
        idx = find(data.burst >= AQDbursts(1) & data.burst <= AQDbursts(2));
        dof = (length(data.P_1)/nfft(i))*2;
        fprintf('ADCP spectra will be averaged with: %0.f Degrees of freedom\n',dof)
    end
    
    %Calculate mean water depth in burst
    [m,n] = size(data.dn(:,idx));
    P = reshape(data.P_1(:,idx),m*n,1)+zp;
    g = zeros(length(P),1);h = zeros(length(P),1);
    for l = 1:length(P)
        g(l,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*P(l,:);
        h(l,:) = ((((-1.82E-15*P(l,:)+2.279E-10)*P(l,:)-2.2512E-5)*P(l,:)+9.72659)*P(l,:))/g(l,:);
    end
    h = nanmean(h);                                 %mean depth during time record
    fprintf('Mean water depth during time record: %0.2f m\n',h)
    
    %Spectral analysis to derive T_i, H_i, and psi_i
    P = zeros(m,n);
    for j = 1:length(idx)
        P(:,j) = detrend(data.P_1(:,idx(j))); %detrend spectrum per segment
    end
    P = reshape(P,m*n,1);                           %redefine P for spectral calcs
    m = nfft(i);                                    %number of windows to average spectra, max(T) output will be m/2
    n = max(size(P));                               %number of data points
    noverlap = round(nfft(i)/4,0);
    
    %%% The following lines from calcspectrum.m %%%
    P = P(:);
    index = 1:m;
    kk = fix((n-noverlap)/(m-noverlap));            %Number of windows
    w = hanning(m);                                 %Window specification
    KMU = kk*norm(w)^2;                             %Normalizing scale factor
    Pxx = zeros(m,1);
    Phx = zeros(m,1);
    %Compute windowed spectra
    for jj=1:kk
        xw = w.*detrend(P(index));
        index = index + (m - noverlap);
        Xx = (abs(fft(xw,nfft(i))).^2)./KMU;
        Ph = angle(fft(xw,nfft(i)));
        Pxx = Pxx + Xx;
        Phx = Phx + Ph;
    end
    % Select first half and let 1st point (DC value) be replaced
    % with 2nd point
    select = [2 2:m/2];
    Pxx = Pxx(select)*2/Fs;
    Phx = Phx(select); %This matches spectra calculated with pwelch
    Cpp = Pxx;
    f = (select - 1)'*Fs/m;
    
    %%% The following lines from waveStats or wds.m %%%
    %Depth correction for pressure sensor elevation
    omega = 2*pi.*f;
    k = qkhf(omega,h)./h;                           %wavenumber from depth
    coshkhz = cosh(k*zp);
    coshkh = cosh(k*h);
    attn = coshkhz./coshkh;                         %transfer function for surface elevation scaling
    attn(attn<0.2) = 0.2;                           %limiter for transfer function
    Spp = Cpp./(attn.^2);                           %surface elevation spectrum
    Hs = 4*sqrt(sum(Spp.*mean(diff(f))));           %Significant wave height
    Tp = 1/(F(Spp == max(Spp(lfc:hfc))));           %peak wave period
    fprintf('Significant wave height during time record: %0.2f m\n',Hs)
    fprintf('Peak wave period during time record: %0.2f s\n\n',Tp)
    pp(i) = plot(f,Spp,'color',cc1(i,:),'linewidth',1.5);
    hold on
    
    if printModelInput == 1
        if contains(data.Gatts.INST_TYPE,'RBR')
            fid = find(f >= 0.03 & f <= 0.75);               %limit output between 33 s and 1.3 s
            f = f(fid);
            Spp = Spp(fid);
            Phx = Phx(fid);
            %%% The following lines adapted from ubspec.m %%%
            amp = sqrt(2*Spp.*mean(diff(f)));               %Wave amp = sqrt(2*S(fi)*df)
            Hi = 2.*amp;                                    %RMS Wave height = 2*amp
            Ti = 1./f;                                      %Wave period per freq
            Psi = wrapTo2Pi(Phx);                           %Wave phase between [0 2pi]
            
            %Print out formatted data blocks to copy/paste into waveProperties
            clc
            fprintf('waveHeights\n%0.0f(\n',length(Hi))
            fprintf('(%0.6f)\n',Hi)
            fprintf(');\n\n')
            
            fprintf('wavePeriods\n%0.0f(\n',length(Ti))
            fprintf('(%0.6f)\n',Ti)
            fprintf(');\n\n')
            
            fprintf('wavePhases\n%0.0f(\n',length(Psi))
            fprintf('(%0.6f)\n',Psi)
            fprintf(');\n\n')
            
            fprintf('waveDirs\n%0.0f(\n',length(Hi))
            fprintf('(%0.0f)\n',zeros(length(Hi),1))
            fprintf(');\n\n')
        end
    end
end
leg = legend(pp,instName,'location','northeast','box','off');
set(gca,'xlim',[0 0.5])
xlabel('\bf\itf (Hz)')
ylabel('\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
prettyfigures('text',11,'labels',12,'box',1,'tlength',[0.005 0.005])

if saveFigs == 1
    export_fig(ff1,[figPath 'MKK18_ModelForcing_' instName{1} '_' instName{2} '_burst' num2str(RBRbursts(1)) 'to' num2str(RBRbursts(2))],'-png','-r600','-nocrop')
    export_fig(ff2,[figPath 'MKK18_AvgSpectra_' instName{1} '_' instName{2} '_burst' num2str(RBRbursts(1)) 'to' num2str(RBRbursts(2))],'-png','-r600','-nocrop')
end