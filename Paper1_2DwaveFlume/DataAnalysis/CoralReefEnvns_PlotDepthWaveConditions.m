%Plot environmental conditions for other USGS reef data to determine
%max wave height, period, and water depths for various reef flats. This
%information supports the modeling effort by determining conditions on
%other reefs worldwide. 
%
% Updates:
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\CheritonFMS2020\';
figPath = 'c:\Users\user\Documents\Data\InstrumentData\Figures\'; %where to save figs
saveFig = 0;
folders = dir(dataPath);
folders = folders(~ismember({folders.name},{'.','..'}));
ff = zeros(length(folders),1);
for i = 1:length(folders) %load data in each data folder
    %Create figure
    ff(i) = figure(i);
    set(ff(i),'PaperOrientation','landscape',...
        'position',[1300 350   1200   500],...
        'color','w','paperpositionmode','auto');
    sp = zeros(3,1);
    fileName = dir([dataPath folders(i).name '\' '*.nc']);
    fileName = {fileName.name};
    cc = brewermap(length(fileName),'Set1');
    pp = zeros(length(fileName),1);
    ps = zeros(length(fileName),1);
    for j = 1:length(fileName)
        disp(['Loading ' fileName{j}])
        data = loadnc([dataPath folders(i).name '\' fileName{j}]);
        
        %work out nfft based on number of samples per burst and a dof of 32
        dof = 32;[spb,~] = size(data.P_1);
        nfft = spb/(dof/2);
        disp('Calculating spectra based on 32 DOF')
        
        %Spectral settings
        lf = 0.05;hf = 0.3; %low and high freq cutoff values [Hz]
        
        %Preallocate variables
        time = zeros(length(data.burst),1);
        Spp = zeros(nfft/2+1,length(data.burst));
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
        
        for jj = 1:length(data.burst)
            
            %Calculate mean water depth in burst
            time(jj) = median(data.dn(:,jj));
            P = data.P_1(:,jj);
            P = cmgbridge(P,100,1000,10000);
            g = zeros(length(P),1);h = zeros(length(P),1);
            for l = 1:length(P)
                g(l,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*P(l,:);
                h(l,:) = ((((-1.82E-15*P(l,:)+2.279E-10)*P(l,:)-2.2512E-5)*P(l,:)+9.72659)*P(l,:))/g(l,:);
            end
            h = mean(h);depth(jj) = h;
            
            %Compute spectra
            P = detrend(P);
            [Cpp,F] = pwelch(P,hanning(nfft),round(nfft*0.7,0),nfft,dt);
            
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
            Spp(:,jj) = Cpp./(attn.^2);                  %surface elevation spectrum
            
            %Calculate Hrms,SS; Hrms, IG
            m0 = sum(Spp(lfc:hfc,jj)*df);                        %zero-th moment
            m1 = sum(F(lfc:hfc).*Spp(lfc:hfc,jj)*df);            %first moment
            m2 = sum((F(lfc:hfc).^2).*Spp(lfc:hfc,jj)*df);       %second moment
            HrmsSS(jj) = 2*sqrt(2.*sum(Spp(lfc:hfc,jj)*df));    %sea-swell rms wave height
            Hs(jj) = 4*sqrt(m0);                               %significant wave height
            Tp(jj) = 1/(F(Spp(:,jj) == max(Spp(lfc:hfc,jj))));   %peak wave period
            Tm(jj) = m0/m1;                                    %mean wave period
            Tz(jj) = sqrt(m0/m2);                              %zero-crossing period
            
            %Infragravity waves
            lfc = find(F >= 0.01,1,'first');            %redefine low freq cutoff
            hfc = find(F <= lf,1,'last');               %redefine high freq cutoff
            HrmsIG(jj) = 2*sqrt(2.*sum(Spp(lfc:hfc,jj)*df));   %infragravity rms wave height
        end
    
        %Plot routine
        sp(1) = subplot(3,1,1);
        plot(time,depth,'color',cc(j,:),'linewidth',1.5)
        hold on
        
        sp(2) = subplot(3,1,2);
        hold on
        pp(j) = plot(time,Hs,'-','color',cc(j,:),'linewidth',1.5);
        
        sp(3) = subplot(3,1,3);
        hold on
        ps(j) = plot(time,smooth(Tp,8),'-','color',cc(j,:),'linewidth',1.5); %smooth is for display purposes
        
        clear data
    end
    instName = regexprep(fileName,'rbr-b.nc','');   
    leg = legend(pp,instName,'location','northeastoutside');

    %Positioning
    set(sp(1),'position',[0.1 0.7 0.75 0.25],...
        'xticklabel',[],... 
        'ylim',[0 8])
    set(sp(2),'position',[0.1 0.4 0.75 0.25],...
        'xticklabel',[],...
        'ylim',[0 4])
    set(sp(3),'position',[0.1 0.1 0.75 0.25],...
        'ylim',[5 25])
    set(leg,'position',[0.87 0.425 0.1 0.2])
    %Labeling
    ylabel(sp(1),'\ith (m)')
    ylabel(sp(2),'\itH_s (m)')
    ylabel(sp(3),'\itT_p (s)')
    datetick(sp(3),'x','mm/dd','keepticks','keeplimits')
    xlabel(sp(3),'Calendar Day')
    title(sp(1),folders(i).name)
end
prettyfigures('text',11,'labels',12,'box',1,'tlength',[0.005 0.005])

