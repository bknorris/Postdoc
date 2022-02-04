%Plot pressure spectra for selected instruments from the 2018 Molokai
%deployment
%
% Updates:
% 01/13/21: Plot wave spectra using raw data files instead of wavestats for
% flexibility
%
% This is version 2 of this script.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'C:\Users\bknorris\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\bknorris\Documents\Data\InstrumentData\Figures\'; %where to save figs
insts = {'MKK18C1601rbr-b.nc';'MKK18HR101aqdHR-b.nc'};
whichRBR = 'C16';
whichAQD = 'HR'; %'HR'; 'LR';
whichBurst = [160 70]; %must be in the same order as insts

ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[1300 350   550   500],...
    'color','w','paperpositionmode','auto');

pp = zeros(length(insts),1);
cc = brewermap(length(insts),'Set1');
for i = 1:length(insts)
    data = loadnc([dataPath insts{i}]);
    
    %Spectral settings
    nfft = [384 768];
    time = zeros(length(data.burst),1);
    Spp = zeros(nfft(i)/2+1,length(data.burst));
    dof = (length(data.P_1)/nfft(i))*2;
    
    %Get attributes from the data
    if contains(data.Gatts.INST_TYPE,'RBR')
        lat = data.Gatts.latitude;                 %latitude for h calcs
        zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
        dt = 1/data.Gatts.sample_interval;         %sample rate [Hz]
        fprintf('RBR spectra will be averaged with: %0.f Degrees of freedom\n',dof)
        
    else
        lat = data.Gatts.latitude;                 %latitude for h calcs
        zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
        dt = data.Gatts.instmeta_AQDSamplingrate;  %sample rate [Hz]
        fprintf('ADCP spectra will be averaged with: %0.f Degrees of freedom\n',dof)
    end
    burstID = find(data.burst == whichBurst(i));
    for j = burstID
        
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
        
        %Compute spectra
        P = detrend(P);
        [Cpp,F] = pwelch(P,hanning(nfft(i)),round(nfft(i)*0.7,0),nfft(i),dt);
        
        %Surface elevation scaling for instrument data
        df = F(3)-F(2);
        omega = 2*pi.*F;
        k = qkhf(omega,h)./h;
        coshkhz = cosh(k*zp);
        coshkh = cosh(k*h);
        attn = coshkhz./coshkh;                     %transfer function for surface elevation scaling
        attn(attn<0.2) = 0.2;                       %limiter for transfer function
        Spp = Cpp./(attn.^2);                       %surface elevation spectrum
    end
    
    pp(i) = plot(F,Spp,'color',cc(i,:),'linewidth',1.5); hold on
end
set(gca,'xlim',[0 0.5])
leg = legend(pp,{'RBR - C16';'ADCP - HR'});
xlabel('\bf\itf (Hz)')
ylabel('\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
title(['\bf\it' sprintf('RBR burst: %0.0f, ADCP burst: %0.0f',whichBurst(1),whichBurst(2))])
prettyfigures('text',11,'labels',12,'box',1)

export_fig(ff1,[figPath 'MKK_Spectra_rbr' whichRBR '_burst' num2str(whichBurst(1)) '_aqd' whichAQD '_burst' num2str(whichBurst(2))],'-png','-r600','-nocrop')


