%Script to analyze the deployment data for the 2018 Molokai deployment
%
% This Script reads in Aquadopp data from either the LR or HR deployment
% above the Molokai Reef (2018 experiment). Burst data are read in, one
% burst at a time, and processed for current and wave statistics. These
% data are saved to disk and optionally plotted. 
%
% This is version 2 of this script.
%
% Changelog: 
% 11/16/2020: added DOF to attributes.
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close;

path = 'C:\Users\bknorris\Documents\Data\InstrumentData\WavesData\';
fname = 'MKK18HR101aqdHR-b.nc'; %'MKK18HR101aqdHR-b.nc'; %'MKK18LR101aqdHR-b.nc';
data = loadnc([path fname]);
ploton = 1; %set to 1 if plots are desired (set to 0 to turn off).

%% Analyze data for current and wave statistics
%Data will be read in burst-by-burst. Output data will be tied to the
%mid-point of each data block to produce a quasi-continuous time series.
disp('Calculating current and wave statistics...')

%Extract metadata from the netCDF file
fs = data.Gatts.instmeta_AQDSamplingrate; %sample rate
lat = data.Gatts.latitude; %latitude -- DOUBLE CHECK THIS!
zp = data.Gatts.transducer_offset_from_bottom; %height of pressure sensor
zuv = zp - data.bindist; %height of velocity bin above bed

%Spectra settings -- CHANGE THESE!
lf = 0.05; %low-freq cutoff
hf = 1.2; %high-freq cutoff
win = 180; %seconds
step = 60; %seconds

%Current and Wave Analysis
avt = fs*step;
nwin = fs*win;
swin = nwin*0.25; %20 second averaging window (to smooth)
DOF = num2str(round((nwin/swin)*2));
disp(['Spectra will be averaged with: ' DOF ' Degrees of freedom'])

wave = struct(); %initialize output data structure
%Populate attribute structure
wave.atts.fileName = fname;
wave.atts.lowFreqCutoff = [num2str(lf) ' Hz'];
wave.atts.hiFreqCutoff = [num2str(hf) ' Hz'];
wave.atts.sampleFreq = [num2str(fs) ' Hz'];
wave.atts.window = [num2str(win) ' seconds'];
wave.atts.step = [num2str(step) ' seconds'];
wave.atts.DOF = DOF;
wave.burst = data.burst;
wave.zuv = zuv;
wave.zp = zp;
tic
for i = 1:length(data.burst) %num bursts
    fprintf('Burst : %0.f\n',data.burst(i))
    p = data.P_1(:,i);
    time = data.dn(:,i);
    %Calculate wave stats in small blocks of length(nwin)
    nsamp = length(time);
    ind = [1 avt:avt:nsamp];
    for j = 1:length(ind)
        if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
            continue
        else
            idx = ind(j):ind(j)+nwin-1;
        end
        %% Derive terms from pressure alone
        t = median(time(idx)); %middle time of data block
        P = p(idx)+zp;
        
        %Bridge short gaps in time series if NaNs are present
        if any(isnan(P))
            disp('Found NaNs in Pressure... interpolating with cmgbridge')
            P = cmgbridge(P,100,100,100);
        end
        %Estimate water depth in each data block
        %SOURCE: Fofonoff and Millard (1982). UNESCO Tech Paper #44.
        g = zeros(length(P),1);h = zeros(length(P),1);
        for l = 1:length(P)
            g(l,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*P(l,:);
            h(l,:) = ((((-1.82E-15*P(l,:)+2.279E-10)*P(l,:)-2.2512E-5)*P(l,:)+9.72659)*P(l,:))/g(l,:);
        end
        h = mean(h); %mean depth during data block
        P = p(idx); %redefine P for spectral calcs
        
        %Detrend time series data for spectral analysis
        P = detrend(P);
        [Cpp,F] = pwelch(P,hanning(nwin),swin*0.7,swin,fs);
        
        %Pressure-derived wave parameters
        df = F(3)-F(2);
        omega = 2*pi.*F;
        k = qkhf(omega,h)./h;
        coshkhz = cosh(k*(h+zp));
        coshkh = cosh(k*h);
        lfc = find(F >= lf,1,'first');              %low freq cutoff
        hfc = find(F <= hf,1,'last');               %high freq cutoff
        attn = coshkhz./coshkh;                     %transfer function for surface elevation scaling
        Spp = Cpp./(attn.^2);                       %surface elevation spectrum
        m0 = sum(Spp(lfc:hfc)*df);                  %zero-th moment
        m1 = sum(F(lfc:hfc).*Spp(lfc:hfc)*df);      %first moment
        m2 = sum((F(lfc:hfc).^2).*Spp(lfc:hfc)*df); %second moment
        Hs = 4*sqrt(m0);                            %significant wave height
        Tp = 1/(F(Spp == max(Spp)));                %peak wave period
        Tm = m0/m1;                                 %mean wave period
        Tz = sqrt(m0/m2);                           %zero-crossing period
        
        %% Save out variables to wave structure
        wave.time(j,i) = t;
        wave.depth(j,i) = h;
        wave.Hs(j,i) = Hs;
        wave.Tp(j,i) = Tp;
        wave.Tm(j,i) = Tm;
        wave.Tz(j,i) = Tz;
        wave.f(:,j,i) = F;
        wave.Cpp(:,j,i) = Cpp; %pressure spectrum at height of instrument (zp)
        wave.Spp(:,j,i) = Spp; %surface pressure spectrum
        
        %% Derive terms from pressure and velocity
        for jj = 1:length(data.depth) %depth bins. NOTE: bottom 2-3 bins are obscured near the bed!
            U = data.u_1205(jj,idx,i)./100; %convert cm/s to m/s;
            V = data.v_1206(jj,idx,i)./100;
            W = data.w_1204(jj,idx,i)./100;
            
            %Bridge short gaps in time series if NaNs are present
            if any([isnan(U) isnan(V) isnan(W)])
                disp('Found NaNs in Velocity... interpolating with cmgbridge')
                fprintf('Bin : %0.f\n',jj)
                U = cmgbridge(U,100,100,100);
                V = cmgbridge(V,100,100,100);
                W = cmgbridge(V,100,100,100);
            end
            
            %Calculate RMS velocities (Luhar et al. 2013)
            %Ec and Nc are mean east and north current velocities
            %Ewrms and Nwrms are rms oscillatory velocities
            %Uc and Uw are total mean oscillatory horizontal velocities
            Ec = (1/nwin)*sum(U);Nc = (1/nwin)*sum(V);
            Ewrms = sqrt((1/nwin)*sum((U-Ec).^2));
            Nwrms = sqrt((1/nwin)*sum((V-Nc).^2));
            Uc = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
            Uw = sqrt(2)*Uwrms;
            
            %Detrend time series data for spectral analysis
            U = detrend(U);
            V = detrend(V);
            W = detrend(W);
            
            %Compute spectra
            [Cuu,~] = pwelch(U,hanning(nwin),swin*0.7,swin,fs);
            [Cvv,~] = pwelch(V,hanning(nwin),swin*0.7,swin,fs);
            [Cpu,~] = cpsd(P,U,hanning(nwin),swin*0.7,swin,fs);
            [Cpv,~] = cpsd(P,V,hanning(nwin),swin*0.7,swin,fs);
            [Cuv,~] = cpsd(U,V,hanning(nwin),swin*0.7,swin,fs);
            Guv = Cuu + Cvv;
            
            %Compute amplitude spectrum for direction/spreading
            Cpu = Cpu.*conj(Cpu);Cpv = Cpv.*conj(Cpv);
            Cuu = Cuu.*conj(Cuu);Cvv = Cvv.*conj(Cvv);
            Cuv = Cuv.*conj(Cuv);
            Dir=57.296*atan2(Cpu(lfc:hfc),Cpv(lfc:hfc));
            Dir=mod(Dir+180,360)-180;                   %Wave direction in [lfc:hfc]; 0 is north
            DirPeak = Dir(Spp == max(Spp));             %Peak wave direction
            R2=((Cuu(lfc:hfc)-Cvv(lfc:hfc)).^2+4*Cuv(lfc:hfc).^2).^.5./(Cuu(lfc:hfc)+Cvv(lfc:hfc));
            Spread = 57.296*((1-R2)/2).^.5;             %Wave spreading
            SpreadPeak = Spread(Spp == max(Spp));       %Peak wave spreading
            omegar = sum(omega(lfc:hfc).*Guv(lfc:hfc)*df)./...
                sum(Guv(lfc:hfc)*df);                   %Orbital wave radian frequency
            ubr = sqrt(2*sum(Guv(lfc:hfc)*df));         %Orbital wave velocity
            Tr = sum(Guv(lfc:hfc).*df)./sum(F(lfc:hfc).*Guv(lfc:hfc).*df); %orbital velocity period
            
            %% Save out variables to wave structure
            wave.Ec(jj,j,i) = Ec;
            wave.Nc(jj,j,i) = Nc;
            wave.Uc(jj,j,i) = Uc;
            wave.Uw(jj,j,i) = Uw;
            wave.dirp(jj,j,i) = DirPeak;
            wave.spreadp(jj,j,i) = SpreadPeak;
            wave.omegar(jj,j,i) = omegar;
            wave.ubr(jj,j,i) = omegar;
            wave.Tr(jj,j,i) = Tr;
        end
    end
end
%% Save data to directory
sfname = regexp(fname,'[^-]+','match');
sfname = [sfname{1} '_waves'];
disp(['Saving ' path sfname '.mat'])
save([path sfname],'wave','-v7.3')
disp(['Analysis completed in: ' num2str(toc/60) ' minutes'])

%% Plot Routine
if ploton == 1
    %Plot the burst pressure data (p)
    F1 = figure(1);
    s(1) = subplot(211);
    title('Pressure')
    ylabel('[DB]')
    hold on
    for i = 1:length(data.burst)
        plot(data.dn(:,i),data.P_1(:,i),'color','b')
    end
    hold off
    %Plot the burst velocity data (u, v, w)
    s(2) = subplot(212);
    title('Velocity')
    ylabel('[cm/s]')
    hold on
    for i = 1:length(data.burst)
        ps(1) = plot(data.dn(:,i),data.u_1205(1,:,i),'color','b');
        ps(2) = plot(data.dn(:,i),data.v_1206(1,:,i),'color','r');
        ps(3) = plot(data.dn(:,i),data.w_1204(1,:,i),'color','g');
    end
    leg = legend(ps,{'U';'V';'W'});
    hold off
    linkaxes(s,'x')
    set(s(1),'xticklabel',[])
    datetickzoom(s(2),'x','mm-dd HH:MM:SS')
    clear s
    
    F2 = figure(2);
    s(1) = subplot(611);
    title('Water Depth')
    ylabel('[m]')
    hold on
    for i = 1:length(wave.burst)
        plot(wave.time(:,i),wave.depth(:,i),'color','b')
        text(median(wave.time(:,i))-0.01,mean(wave.depth(:,i))+0.05,['\bf\it' num2str(wave.burst(i))])
    end
    hold off
    s(2) = subplot(612);
    title('Significant Wave Height')
    ylabel('[m]')
    hold on
    for i = 1:length(wave.burst)
        plot(wave.time(:,i),wave.Hs(:,i),'color','b')
    end
    hold off
    s(3) = subplot(613);
    title('Mean Wave Period')
    ylabel('[s]')
    hold on
    for i = 1:length(wave.burst)
        plot(wave.time(:,i),wave.Tm(:,i),'color','b')
    end
    hold off
    
    s(4) = subplot(614);
    title('Mean Current Velocity')
    ylabel('[m/s]')
    cc = brewermap(15,'Blues');
    hold on
    for i = 1:length(wave.burst)
        for j = 1:15
            plot(wave.time(:,i),wave.Uc(j,:,i),'color',cc(j,:))
        end
    end
    hold off
    s(5) = subplot(615);
    title('Mean Wave Velocity')
    ylabel('[m/s]')
    hold on
    for i = 1:length(wave.burst)
        for j = 1:15
            plot(wave.time(:,i),wave.Uw(j,:,i),'color',cc(j,:))
        end
    end
    hold off
    s(6) = subplot(616);
    title('Peak Wave Direction')
    ylabel('[deg]')
    hold on
    for i = 1:length(wave.burst)
        for j = 1:15
            plot(wave.time(:,i),wave.dirp(j,:,i),'color',cc(j,:))
        end
    end
    hold off
    set(s(6),'ylim',[0 360],'ytick',0:90:360)
    
    linkaxes(s,'x')
    set([s(1) s(2) s(3) s(4) s(5)],'xticklabel',[])
    datetickzoom(s(6),'x','mm-dd HH:MM:SS')
    prettyfigures('text',12,'labels',13,'box',1,'tickdir','in','tlength',[0.005 0.005])
end