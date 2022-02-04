%Plot environmental conditions for the RBR deployments to determine min and
%max wave height, period, and water depths for modeling.
% Updates:
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc
dataPath = 'C:\Users\user\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\user\Documents\Data\InstrumentData\Figures\'; %where to save figs
insts = {'MKK18C1501rbr-b.nc';'MKK18C1601rbr-b.nc'};
instName = {'C15';'C16'};
figName = 'MKK18_DepthWaveConditions_C15_C16';
startTime = datenum(2018,06,25,00,00,00);
endTime = datenum(2018,08,10,00,00,00);
saveFig = 0;

%Create figure
ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[1300 350   700   600],...
    'color','w','paperpositionmode','auto');
sp = zeros(5,1);
pp = zeros(2,1);
ps = zeros(2,1);
cc = [34 94 168;29 145 192]./255;
for i = 1:length(insts)
    %Load data
    data = loadnc([dataPath insts{i}]);

    %Spectral settings
    nfft = [256 256];
    dof = (length(data.P_1)/nfft(i))*2;
    lf = 0.05;hf = 0.3; %low and high freq cutoff values [Hz]
    
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
        h = mean(h);depth(j) = h;
        
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
        HrmsSS(j) = 2*sqrt(2.*sum(Spp(lfc:hfc,j)*df));    %sea-swell rms wave height
        Hs(j) = 4*sqrt(m0);                               %significant wave height
        Tp(j) = 1/(F(Spp(:,j) == max(Spp(lfc:hfc,j))));   %peak wave period
        Tm(j) = m0/m1;                                    %mean wave period
        Tz(j) = sqrt(m0/m2);                              %zero-crossing period
        
        %Infragravity waves
        lfc = find(F >= 0.01,1,'first');            %redefine low freq cutoff
        hfc = find(F <= lf,1,'last');               %redefine high freq cutoff
        HrmsIG(j) = 2*sqrt(2.*sum(Spp(lfc:hfc,j)*df));   %infragravity rms wave height  
    end
    
    %Plot routine
    sp(1) = subplot(3,1,1);
    plot(time,detrend(depth),'color',cc(i,:),'linewidth',1.5)
    hold on
    
    sp(2) = subplot(3,1,2);
    hold on
    pp(i) = plot(time,Hs,'-','color',cc(i,:),'linewidth',1.5);
    
    sp(3) = subplot(3,1,3);
    hold on
    ps(i) = plot(time,smooth(Tp,8),'-','color',cc(i,:),'linewidth',1.5); %smooth is for display purposes
    
    if contains(insts{i},'C16')
        fprintf('Min water depth at C16: %0.2f m\n',min(min(depth)))
        fprintf('Max water depth at C16: %0.2f m\n',max(max(depth)))
        fprintf('Min wave height at C16: %0.2f m\n',min(min(Hs)))
        fprintf('Min wave height at C16: %0.2f m\n',max(max(Hs)))
        fprintf('Min wave period at C16: %0.2f m\n',min(min(smooth(Tp,8))))
        fprintf('Max wave period at C16: %0.2f m\n',max(max(smooth(Tp,8))))
    end
end

%Positioning
set(sp(1),'position',[0.12 0.7 0.75 0.25],...
    'xticklabel',[],'xlim',[startTime endTime],...
    'xtick',startTime:datenum(0,0,5,0,0,0):endTime,...
    'ylim',[-0.5 0.5])
set(sp(2),'position',[0.12 0.4 0.75 0.25],...
    'xticklabel',[],'xlim',[startTime endTime],...
    'xtick',startTime:datenum(0,0,5,0,0,0):endTime,...
    'ylim',[0 0.8])
set(sp(3),'position',[0.12 0.1 0.75 0.25],...
    'xlim',[startTime endTime],...
    'xtick',startTime:datenum(0,0,5,0,0,0):endTime,...
    'ylim',[5 25])
leg1 = legend(pp,instName,'location','northeast');
leg2 = legend(ps,instName,'location','northeast');
%Labeling
ylabel(sp(1),'\ith (m)')
ylabel(sp(2),'\itH_s (m)')
ylabel(sp(3),'\itT_p (s)')
datetick(sp(3),'x','mm/dd','keepticks','keeplimits')
xlabel(sp(3),'Calendar Day')
prettyfigures('text',11,'labels',12,'box',1,'tlength',[0.005 0.005])
set(leg1,'position',[0.78 0.59 0.05 0.05],'box','off')
set(leg2,'position',[0.78 0.29 0.05 0.05],'box','off')
if saveFig == 1
    export_fig(ff1,[figPath  'MKK18_RepresentativeConditions_C15_C16_V2'],'-png','-r600','-nocrop')
end

