%Plot pressure spectra for selected instruments from the 2018 Molokai
%to visualized  representative periods of time.
%
% Updates:
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'C:\Users\user\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\user\Documents\Data\InstrumentData\Figures\'; %where to save figs
insts = {'MKK18C1501rbr-b.nc';'MKK18C1601rbr-b.nc'};
instName = {'C15';'C16'};
figName = 'MKK18_RepresentativeSpectra_C15_C16';
startTime = datenum(2018,06,25,00,00,00);
endTime = datenum(2018,07,25,00,00,00);
saveFig = 0;

%Create figure
ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[1300 350   1050   800],...
    'color','w','paperpositionmode','auto');
sp = zeros(5,1);
pp = zeros(2,1);
cc1 = [34 94 168;29 145 192]./255;
cc2 = brewermap(300,'*YlGnBu');
for i = 1:length(insts)
    %Load data
    data = loadnc([dataPath insts{i}]);

    %Spectral settings
    nfft = [256 256];
    dof = (length(data.P_1)/nfft(i))*2;
    lf = 0.03;hf = 0.3; %low and high freq cutoff values [Hz]
    
    %Preallocate variables
    time = zeros(length(data.burst),1);
    Spp = zeros(nfft(i)/2+1,length(data.burst));
    depth = zeros(length(data.burst),1);
    Hs = zeros(length(data.burst),1);
    HrmsSS = zeros(length(data.burst),1);
    HrmsIG = zeros(length(data.burst),1);
    Tp = zeros(length(data.burst),1);
    Tm = zeros(length(data.burst),1);
    Tz = zeros(length(data.burst),1);
    
    %Get attributes from the data
    lat = data.Gatts.latitude;                 %latitude for h calcs
    zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
    dt = 1/data.Gatts.sample_interval;         %sample rate [Hz]

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
        h = mean(h);depth(:,j) = h;
        
        %Compute spectra
        P = detrend(P);
        [Cpp,F] = pwelch(P,hanning(nfft(i)),round(nfft(i)*0.7,0),nfft(i),dt);
        
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
        
        %Calculate Hrms,SS; Hrms, IG
        m0 = sum(Spp(lfc:hfc,j)*df);                        %zero-th moment
        m1 = sum(F(lfc:hfc).*Spp(lfc:hfc,j)*df);            %first moment
        m2 = sum((F(lfc:hfc).^2).*Spp(lfc:hfc,j)*df);       %second moment
        HrmsSS(:,j) = 2*sqrt(2.*sum(Spp(lfc:hfc,j)*df));    %sea-swell rms wave height
        Hs(:,j) = 4*sqrt(m0);                               %significant wave height
        Tp(:,j) = 1/(F(Spp(:,j) == max(Spp(lfc:hfc,j))));   %peak wave period
        Tm(:,j) = m0/m1;                                    %mean wave period
        Tz(:,j) = sqrt(m0/m2);                              %zero-crossing period
        
        %Infragravity waves
        lfc = find(F >= 0.01,1,'first');            %redefine low freq cutoff
        hfc = find(F <= lf,1,'last');               %redefine high freq cutoff
        HrmsIG(:,j) = 2*sqrt(2.*sum(Spp(lfc:hfc,j)*df));   %infragravity rms wave height  
    end
    
    %Plot routine
    sp(1) = subplot(5,1,1);
    pp(i) = plot(time(:,:),depth(i,:),'color',cc1(i,:),'linewidth',1.5);
    hold on
    
    sp(2) = subplot(5,1,2);
    hold on
    ps(1) = plot(time(:,:),Hs(i,:),'-','color',cc1(i,:),'linewidth',1.5);
    ps(2) = plot(time(:,:),HrmsSS(i,:),'+','color',cc1(i,:),'linewidth',1.5);
    ps(3) = plot(time(:,:),HrmsIG(i,:),'.','color',cc1(i,:),'linewidth',1.5);
    if i == 2
        legend(ps,{'H_s';'H_{rms,SS}';'H_{rms,IG}'},'location','eastoutside')
    end
    clear ps
    
    sp(3) = subplot(5,1,3);
    hold on
    ps(1) = plot(time(:,:),smooth(Tp(i,:),5),'-','color',cc1(i,:),'linewidth',1.5); %smooth is for display purposes
    ps(2) = plot(time(:,:),Tm(i,:),'+','color',cc1(i,:),'linewidth',1.5);
    ps(3) = plot(time(:,:),Tz(i,:),'.','color',cc1(i,:),'linewidth',1.5);
    if i == 2
        legend(ps,{'T_p';'T_m';'T_z'},'location','eastoutside')
    end
    
    sp(i+3) = subplot(5,1,i+3);
    imagesc(time(:,:),F,squeeze(Spp(:,:)))
    caxis([0.01 0.05])
    colormap(cc2);
    title(['\bf\it' instName{i}])
end
legend(pp,instName,'location','eastoutside')

%Positioning
set(sp(1),'position',[0.12 0.85 0.75 0.12],...
    'xticklabel',[],'xlim',[startTime endTime],...
    'xtick',startTime:datenum(0,0,1,0,0,0):endTime,...
    'ylim',[0 2])
set(sp(2),'position',[0.12 0.69 0.75 0.12],...
    'xticklabel',[],'xlim',[startTime endTime],...
    'xtick',startTime:datenum(0,0,1,0,0,0):endTime,...
    'ylim',[0 0.5])
set(sp(3),'position',[0.12 0.52 0.75 0.12],...
    'xticklabel',[],'xlim',[startTime endTime],...
    'xtick',startTime:datenum(0,0,1,0,0,0):endTime,...
    'ylim',[5 20])
set(sp(4),'xlim',[startTime endTime],...
    'xticklabel',[],...
    'xtick',startTime:datenum(0,0,1,0,0,0):endTime,...
    'ylim',[0 0.5],...
    'position',[0.12 0.3 0.75 0.18])
set(sp(5),'xlim',[startTime endTime],...
    'xtick',startTime:datenum(0,0,1,0,0,0):endTime,...
    'ylim',[0 0.5],...
    'position',[0.12 0.08 0.75 0.18])

%Labeling
ylabel(sp(1),'\bf\itDepth (m)')
ylabel(sp(2),'\bf\itWave Height (m)')
ylabel(sp(3),'\bf\itWave Period (s)')
ylabel(sp(4),'\bf\itf (Hz)')
ylabel(sp(5),'\bf\itf (Hz)')
datetick(sp(5),'x','mm-dd','keepticks','keeplimits')
xlabel(sp(5),'\bf\itDate in 2018')
cb = colorbar;
ylabel(cb,'\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
set(cb,'position',[0.89 0.08 0.015 0.4],'linewidth',1.5,'tickdir','out')
prettyfigures('text',11,'labels',12,'box',1,'tlength',[0.005 0.005])
if saveFig == 1
    export_fig(ff1,[figPath  'MKK18_RepresentativeConditions_C15_C16_V2'],'-png','-r600','-nocrop')
end

