%Script to analyze the deployment data for the 2018 Molokai deployment
%
% This Script reads in RBR data from the C-line of the MKK deployment one
% at a time and analyzes the pressure time series for wave statistics (Hs,
% Tp, Tm, Tz, and pressure spectra), then writes out a 'wave' file which
% contains these derived quantities.
%
% This is version 2 of this script.
%
% Changelog:
% 11/16/2020: added DOF to attributes, rms wave heights for SS and IG bands
% 11/17/2020: re-wrote script to do burst processing more efficiently
%
% BKN - USGS PCMSC 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% close all;

path = 'C:\Users\bknorris\Documents\Data\InstrumentData\WavesData\';
files = dir([path '*rbr*' '.nc']);
fileList = {files.name};

%Aquadopp start and end times will be written to the data file for later
%referece
startTime = datenum('22-Jun-2018 23:01:29');
endTime = datenum('25-Jun-2018 18:32:29');

%Spectra settings -- CHANGE THESE!
lf = 0.05; %low-freq cutoff
hf = 1.2; %high-freq cutoff

%% Analyze data for current and wave statistics
%Data will be read in burst-by-burst. Output data will be tied to the
%mid-point of each data block to produce a quasi-continuous time series.
t1 = tic;
for i = 1:length(fileList)
    disp(['Loading ' fileList{i}])
    data = loadnc([path fileList{i}]);
    t2 = tic;
    disp('Calculating current and wave statistics...')
    
    %Extract metadata from the netCDF file
    fs = 1/data.Gatts.sample_interval; %sample rate = 1/sample interval
    lat = data.Gatts.latitude; %latitude
    zp = data.Gatts.initial_instrument_height; %height of pressure sensor
    spb = data.Gatts.samples_per_burst;
    mode = data.Gatts.mode;
    
    %If continuous, subdivide t-s into "bursts"
    if strcmp(mode,'continuous')
        nwin = round(spb/10,0);
    elseif strcmp(mode,'burst')
        nwin = spb;
    end
    swin = nwin.*0.25;
    DOF = num2str(round((nwin/swin)*2));
    disp(['Spectra will be averaged with: ' DOF ' Degrees of freedom'])
    
    %Populate attribute structure
    wave = struct(); %initialize output data structure
    wave.atts.fileName = fileList{i};
    wave.atts.lowFreqCutoff = [num2str(lf) ' Hz'];
    wave.atts.hiFreqCutoff = [num2str(hf) ' Hz'];
    wave.atts.sampleFreq = [num2str(fs) ' Hz'];
    wave.atts.window = [num2str(nwin/fs) ' seconds'];
    wave.atts.DOF = str2double(DOF);
    wave.atts.AQDPstart = startTime;
    wave.atts.AQDPend = endTime;
    wave.zp = zp;
    fprintf('Processing %0.f bursts\n',length(data.burst))
    for j = 1:length(data.burst) %num bursts
        wave.burst(j) = data.burst(j);
        p = data.P_1(:,j);
        time = data.dn(:,j);
        nsamp = length(time);
        idx = [1 nwin:nwin:nsamp];
        for jj = 1:length(idx)-1
            t = median(time(idx(jj):idx(jj+1))); %middle time of data block
            P = p(idx(jj):idx(jj+1))+zp;
            
            %Estimate water depth in each data block
            %SOURCE: Fofonoff and Millard (1982). UNESCO Tech Paper #44.
            g = zeros(length(P),1);h = zeros(length(P),1);
            for l = 1:length(P)
                g(l,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*P(l,:);
                h(l,:) = ((((-1.82E-15*P(l,:)+2.279E-10)*P(l,:)-2.2512E-5)*P(l,:)+9.72659)*P(l,:))/g(l,:);
            end
            h = nanmean(h); %mean depth during data block
            P = p(idx(jj):idx(jj+1)); %redefine P for spectral calcs
            
            %Bridge short gaps in time series if NaNs are present
            if any(isnan(P))
                fprintf('Burst : %0.f\n',data.burst(j))
                disp('Found NaNs in Pressure... interpolating with cmgbridge')
                P = cmgbridge(P,100,100,1000);
            end
            
            %Detrend time series data for spectral analysis
            P = detrend(P);
            [Cpp,F,Cppx] = pwelch(P,hanning(nwin),swin*0.7,swin,fs);
            
            %Pressure-derived wave parameters
            df = F(3)-F(2);
            omega = 2*pi.*F;
            k = qkhf(omega,h)./h;
            coshkhz = cosh(k*zp);
            coshkh = cosh(k*h);
            lfc = find(F >= lf,1,'first');              %low freq cutoff
            hfc = find(F <= hf,1,'last');               %high freq cutoff
            attn = coshkhz./coshkh;                     %transfer function for surface elevation scaling
            attn(attn<0.2) = 0.2;                       %limiter for transfer function
            Spp = Cpp./(attn.^2);                       %surface elevation spectrum
            m0 = sum(Spp(lfc:hfc)*df);                  %zero-th moment
            m1 = sum(F(lfc:hfc).*Spp(lfc:hfc)*df);      %first moment
            m2 = sum((F(lfc:hfc).^2).*Spp(lfc:hfc)*df); %second moment
            HrmsSS = 2*sqrt(2.*sum(Spp(lfc:hfc)*df));   %sea-swell rms wave height
            Hs = 4*sqrt(m0);                            %significant wave height
            Tp = 1/(F(Spp == max(Spp(lfc:hfc))));       %peak wave period
            Tm = m0/m1;                                 %mean wave period
            Tz = sqrt(m0/m2);                           %zero-crossing period
            
            %Infragravity waves
            lfc = find(F >= 0.01,1,'first');            %redefine low freq cutoff
            hfc = find(F <= lf,1,'last');               %redefine high freq cutoff
            HrmsIG = 2*sqrt(2.*sum(Spp(lfc:hfc)*df));   %infragravity rms wave height
            
            %% Save out variables to wave structure
            wave.time(jj,j) = t;
            wave.depth(jj,j) = h;
            wave.HrmsSS(jj,j) = HrmsSS;
            wave.HrmsIG(jj,j) = HrmsIG;
            wave.Hs(jj,j) = Hs;
            wave.Tp(jj,j) = Tp;
            wave.Tm(jj,j) = Tm;
            wave.Tz(jj,j) = Tz;
            wave.f(:,jj,j) = F;
            wave.Cpp(:,jj,j) = Cpp; %pressure spectrum at height of instrument (zp)
            wave.Spp(:,jj,j) = Spp; %surface pressure spectrum
        end
    end
    %% Save data to directory
    sfname = regexp(fileList{i},'[^-]+','match');
    sfname = [sfname{1} '_waves'];
    disp(['Saving ' path sfname '.mat'])
    save([path sfname],'wave','-v7.3')
    disp(['File processed in: ' num2str(toc(t2)/60) ' minutes'])
end
disp(['Analysis completed in: ' num2str(toc(t1)/60) ' minutes'])
