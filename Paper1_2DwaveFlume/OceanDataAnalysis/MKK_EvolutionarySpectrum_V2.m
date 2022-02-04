%Compare pressure spectra from C15 and the ADCPs for the entire time
%series. This script plots an evolutionary spectrum from both instruments
%for direct comparison.
%
%
% This is version 2 of this script.
%
% Changelog:
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all
dataPath = 'c:\Users\user\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\user\Documents\Data\InstrumentData\Figures\'; %where to save figs
insts = {'MKK18C1501rbr-b.nc';'MKK18HR101aqdHR-b.nc'};
whichAQD = 'HR'; %'HR'; 'LR';

ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[1300 350   950   500],...
    'color','w','paperpositionmode','auto');

sp = zeros(length(insts),1);
m = zeros(1,1);
startTime = zeros(length(insts),1);
endTime = zeros(length(insts),1);
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
        RBRtime = zeros(1,length(data.burst));     %save out times for analysis
        RBRfreq = zeros(nfft(i)/2+1,length(data.burst)); %save out frequencies
        RBRspec = zeros(nfft(i)/2+1,length(data.burst)); %save out surface spectrum
        RBRburst = zeros(1,length(data.burst));          %save out burst no for plotting
        fprintf('RBR spectra will be averaged with: %0.f Degrees of freedom\n',dof)

    else
        lat = data.Gatts.latitude;                 %latitude for h calcs
        zp = data.Gatts.initial_instrument_height; %height of pressure sensor [m]
        dt = data.Gatts.instmeta_AQDSamplingrate;  %sample rate [Hz]
        AQDtime = zeros(1,length(data.burst));     %save out times for analysis
        AQDfreq = zeros(nfft(i)/2+1,length(data.burst)); %save out frequencies
        AQDspec = zeros(nfft(i)/2+1,length(data.burst)); %save out surface spectrum
        AQDburst = zeros(1,length(data.burst));          %save out burst no for plotting
        fprintf('ADCP spectra will be averaged with: %0.f Degrees of freedom\n',dof)
    end
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
        Spp(:,j) = Cpp./(attn.^2);                  %surface elevation spectrum
        
        if contains(data.Gatts.INST_TYPE,'RBR')
            RBRtime(j) = time(j);
            RBRfreq(:,j) = F;
            RBRspec(:,j) = Spp(:,j);
            RBRburst(j) = data.burst(j);
        else
            AQDtime(j) = time(j);
            AQDfreq(:,j) = F;
            AQDspec(:,j) = Spp(:,j);
            AQDburst(j) = data.burst(j);
        end
    end
    startTime(i) = min(time);
    endTime(i) = max(time);
    
    %Plot routine
    sp(i) = subplot(1,2,i);
    imagesc(time(:,:),F,squeeze(Spp(:,:)))
    datetickzoom('x','mm-dd','keepticks','keeplimits')
    caxis([0.01 0.05])
    title(['\bf\it' insts{i}])
    hold on
end
ylabel(sp(1),'\bf\itf (Hz)')
xlabel(sp(1),'\bf\itDate in 2018')
xlabel(sp(2),'\bf\itDate in 2018')
cb = colorbar;
ylabel(cb,'\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
set(sp(1),'xlim',[startTime(1) datenum(2018,06,29,0,0,0)],...
    'xtick',startTime(1):datenum(0,0,1,0,0,0):datenum(2018,06,29,0,0,0),...
    'ylim',[0 0.5],...
    'position',[0.1 0.12 0.45 0.8])
set(sp(2),'xlim',[startTime(2) endTime(2)],...
    'ylim',[0 0.5],...
    'position',[0.59 0.12 0.25 0.8],...
    'yticklabel',[],'xcolor','r','ycolor','r')
set(cb,'position',[0.86 0.12 0.025 0.8],'linewidth',1.5,'tickdir','out')
prettyfigures('text',11,'labels',12,'box',1)

%Plot rectangle on RBR subplot to denote ADCP times
rx = [startTime(2) endTime(2) endTime(2) startTime(2) startTime(2)];
ry = [0 0 0.5 0.5 0];
plot(sp(1),rx,ry,'-r','linewidth',1.5)

%Select points on the ADCP record for plotting
yes = 0;
repeat = 1;
while repeat == 1
    while yes ~= 1
        disp('USER: select a time in the ADCP record to plot spectrum')
        jj = waitforbuttonpress;
        point1 = get(gca,'CurrentPoint');
        %Plot selected point as a green line on the figure
        line1 = plot(sp(2),point1(1)*ones(10,1),linspace(0,0.5,10),'-g','linewidth',1.5);
        %Locate selected time in RBR time series
        timeDiff = etime(datevec(RBRtime(1)),datevec(point1(1)));
        if timeDiff > 0
            disp('USER selected a time that is outside of the RBR''s deployment window. Pick another.')
            yes = 0;
            delete(line1);
        elseif timeDiff < 0
            disp('Plotting Spectra...')
            [~,RBidx] = min(abs(RBRtime-point1(1)));
            [~,ADidx] = min(abs(AQDtime-point1(1)));
            fprintf('Burst time: %s\n',datestr(AQDtime(ADidx)))
            line2 = plot(sp(1),RBRtime(RBidx)*ones(10,1),linspace(0,0.5,10),'-g','linewidth',1.5);
            ff2 = figure(2);
            set(ff2,'PaperOrientation','landscape',...
                'position',[700 350   550   500],...
                'color','w','paperpositionmode','auto');
            pp(1) = plot(RBRfreq(:,RBidx),RBRspec(:,RBidx),'-b','linewidth',1.5); hold on
            pp(2) = plot(AQDfreq(:,ADidx),AQDspec(:,ADidx),'-r','linewidth',1.5);
            set(gca,'xlim',[0 0.5])
            leg = legend(pp,{'RBR';'ADCP'});
            xlabel('\bf\itf (Hz)')
            ylabel('\bf\itS(f)_{\eta\eta}  [m^2/Hz]')
            title(['\bf\it' sprintf('RBR burst: %0.0f, ADCP burst: %0.0f',RBRburst(RBidx),AQDburst(ADidx))])
            prettyfigures('text',11,'labels',12,'box',1)
            prompt = 'Is this selection correct? [y/n] ';
            str = input(prompt,'s');
            if strcmp(str,'y')
                yes = 1;
                repeat = 0;
            elseif strcmp(str,'n')
                yes = 0;
                delete(line1);
                delete(line2);
                close(ff2);
            end
        end
    end
end
disp('Saving figures...')
export_fig(ff1,[figPath 'MKK_EvoSpectra_rbr_burst' num2str(RBRburst(RBidx)) '_aqd' whichAQD '_burst' num2str(AQDburst(ADidx))],'-png','-r600','-nocrop')
export_fig(ff2,[figPath 'MKK_Spectra_rbr_burst' num2str(RBRburst(RBidx)) '_aqd' whichAQD '_burst' num2str(AQDburst(ADidx))],'-png','-r600','-nocrop')

