%Compare pressure spectra from C15 and the ADCPs for the entire time
%series. This script plots an evolutionary spectrum from both instruments
%for direct comparison.
%
%
% This is version 1 of this script.
%
% Changelog:
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
insts = {'MKK18C1501rbr-b.nc';'MKK18HR101aqdHR-b.nc'};

ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[500 350   950   500]);
sp = zeros(length(insts),1);
m = zeros(1,1);
for i = 1:length(insts)
    data = loadnc([dataPath insts{i}]);
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
    %Spectral settings
    nfft = [256 512];
    time = zeros(length(data.burst),1);
    Spp = zeros(nfft(i)/2+1,length(data.burst));
    
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
        h = mean(h);
        
        %% Compute spectra
        P = detrend(P);
        [Cpp,F] = pwelch(P,hanning(nfft(i)),round(nfft(i)*0.7,0),nfft(i),dt);
        
        %Do surface elevation scaling for the instrument data
        df = F(3)-F(2);
        omega = 2*pi.*F;
        k = qkhf(omega,h)./h;
        coshkhz = cosh(k*zp);
        coshkh = cosh(k*h);
        attn = coshkhz./coshkh;                     %transfer function for surface elevation scaling
        attn(attn<0.2) = 0.2;                       %limiter for transfer function
        Spp(:,j) = Cpp./(attn.^2);                       %surface elevation spectrum
    end
    sp(i) = subplot(1,2,i);
    imagesc(time,F,Spp)
    datetickzoom('x','mm-dd','keepticks','keeplimits')
    caxis([0.01 0.05])
    set(gca,'xlim',[min(time) max(time)],'ylim',[0 0.5])
    title(['\bf\it' insts{i}])
end
ylabel(sp(1),'\bf\itf (Hz)')
xlabel(sp(1),'\bf\itDate in 2018')
xlabel(sp(2),'\bf\itDate in 2018')
cb = colorbar;
ylabel(cb,'\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
prettyfigures('text',11,'labels',12,'box',1)

