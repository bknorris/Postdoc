%Calculate wave energy flux between the C15 RBR and the HR Aquadopp to
%determine the total wave dissipation due to breaking (epsilon_b) and
%bottom friction (epsilon_f). Wave energy flux calculations are based off
%of Huang et al., (2012) and Lowe et al., 2005. The wave dissipation rate is
%
%D = deltaF/(deltaX*cos(theta)),
%
%where deltaF is the change in energy flux between two measurement points,
%deltaX is the distance between measurement points, and cos(theta) is the
%direction of wave propagation.
%
% Updates:
%
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
dataPath = 'c:\Users\bknorris\Documents\Data\InstrumentData\WavesData\';
figPath = 'c:\Users\bknorris\Documents\Data\InstrumentData\Figures\'; %where to save figs
insts = {'MKK18C1501rbr-b.nc';'MKK18HR101aqdHR-b.nc'};
instName = {'C15';'HR-Aquadopp'};
RBRbursts = [151 158]; %RBR bursts to analyze (range: 141 - 162)
AQDbursts = [61 68]; %ADCP bursts to analyze (range: 51 - 72)
nfft = [128 128]; %spectral window settings for instruments
win = 120; %5 minute averaging interval for wave energy flux
saveFigs = 0; %save figures; 1 = yes, 0 = no
%% Process instrument data
if length(RBRbursts) == 1
    burstTxt = sprintf('%d',RBRbursts(1));
else
    burstTxt = sprintf('%d-%d',RBRbursts(1),RBRbursts(2));
end

%Load both instruments
RBR = loadnc([dataPath insts{1}]);
AQD = loadnc([dataPath insts{2}]);

%Extract important values from RBR file first
lat = RBR.Gatts.latitude;                 %latitude for h calcs
RBRzp = RBR.Gatts.initial_instrument_height; %height of pressure sensor [m]
fs = 1/RBR.Gatts.sample_interval;         %sample interval [s]
idx1 = find(RBR.burst >= RBRbursts(1) & RBR.burst <= RBRbursts(2));
dof = (length(RBR.P_1)/nfft(1))*2;
fprintf('RBR spectra will be averaged with: %0.f Degrees of freedom\n',dof)

[m,n] = size(RBR.dn(:,idx1));
P1 = reshape(RBR.P_1(:,idx1),m*n,1);
RBRtime = reshape(RBR.dn(:,idx1),m*n,1);

%Now do the ADCP
AQDzp = AQD.Gatts.initial_instrument_height; %height of pressure sensor [m]
idx1 = find(AQD.burst >= AQDbursts(1) & AQD.burst <= AQDbursts(2));

[m,n] = size(AQD.dn(:,idx1));
P2 = reshape(AQD.P_1(:,idx1),m*n,1);
U = reshape(AQD.u_1205(1,:,idx1),m*n,1)./100; %convert from cm/s to m/s; using the top bin
V = reshape(AQD.v_1206(1,:,idx1),m*n,1)./100; %convert from cm/s to m/s; using the top bin
AQDtime = reshape(AQD.dn(:,idx1),m*n,1);

%Sync measurments between the RBR and Aquadopp
if (RBRtime(1)-AQDtime(1)) == 0 %check this
    disp('Start times for the RBR and Aquadopp are the same')
else
    error('Need to ensure RBR and Aquadopp start times are the same!')
end
P2 = interp1(AQDtime,P2,RBRtime); %interp Aquadopp pressure to same length as the RBR time series
U = interp1(AQDtime,U,RBRtime);
V = interp1(AQDtime,V,RBRtime);

%Put NaNs in where the gaps in bursts are
RBRnumBursts = RBRbursts(2)-RBRbursts(1);
RBRburstLength = length(RBR.dn);
RBRburstNaNs = RBRburstLength+1:RBRburstLength:RBRnumBursts*RBRburstLength+1;
RBRtime(RBRburstNaNs) = NaN;

%% Compute WEF
avt = fs*win;
idx2 = 1:avt:length(P1);
%Initialize variables
F_rbr = NaN(length(idx2)-1,1);
F_aqd = NaN(length(idx2)-1,1);
theta_aqd = NaN(length(idx2)-1,1);
time = NaN(length(idx2)-1,1);
for i = 1:length(idx2)-1
    %Do the RBR first
    inds = idx2(i):idx2(i+1);
    pp = P1(inds);
    tt = RBRtime(inds);
    if any(find(isnan(tt))) == 1
        continue
    else
        rho = 1023; %density of seawater at 25C (kg/m^3)
        g = 9.81; %gravitational constant (m/s^2)
        
        h = zeros(length(pp),1);
        for l = 1:length(pp)
            gg = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*pp(l,:);
            h(l,:) = ((((-1.82E-15*pp(l,:)+2.279E-10)*pp(l,:)-2.2512E-5)*pp(l,:)+9.72659)*pp(l,:))/gg;
        end
        h = nanmean(h);                             %mean depth during time record
        [Cpp,F] = pwelch(detrend(pp),hanning(nfft(1)),round(nfft(1)*0.25,0),nfft(1),fs);
        
        %Compute wave energy flux
        zp = RBR.Gatts.initial_instrument_height; %height of pressure sensor [m]
        df = F(3)-F(2);
        omega = 2*pi.*F;
        k = qkhf(omega,h)./h;
        kh = k*h;
        coshkhz = cosh(k*zp);
        coshkh = cosh(k*h);
        attn = coshkhz./coshkh;                     %transfer function for surface elevation scaling
        attn(attn<0.2) = 0.2;                       %limiter for transfer function
        Spp = Cpp./(attn.^2);                       %surface elevation spectrum
        amp = sqrt(2*Spp.*df);                      %Wave amp = sqrt(2*S(fi)*df)
        E = 0.5*rho*g*amp.^2;                       %Wave energy [Lowe et al., 2005]
        Cg = 0.5*(1+((2*kh)./(sinh(2*kh)))).*(omega./k); %Wave group velocity [Lowe et al., 2005]
        F_rbr(i) = trapz(E(2:end).*Cg(2:end));      %Total Wave energy flux (W/m)
        time(i) = tt(180);                          %Datetime for windowed interval
        
        %Now do the Aquadopp
        pp = P2(inds);
        uu = U(inds);
        vv = V(inds);
        rho = 1023; %density of seawater at 25C (kg/m^3)
        g = 9.81; %gravitational constant (m/s^2)
        
        h = zeros(length(pp),1);
        for l = 1:length(pp)
            gg = 9.780318*(1.0+(5.2788E-3+2.36E-5*lat)*lat)+1.092E-6*pp(l,:);
            h(l,:) = ((((-1.82E-15*pp(l,:)+2.279E-10)*pp(l,:)-2.2512E-5)*pp(l,:)+9.72659)*pp(l,:))/gg;
        end
        h = nanmean(h);                             %mean depth during time record
        [Cpp,F] = pwelch(detrend(pp),hanning(nfft(2)),round(nfft(2)*0.25,0),nfft(2),fs);
        [Cpu,~] = cpsd(detrend(pp),detrend(uu),hanning(nfft(2)),round(nfft(2)*0.25,0),nfft(2),fs);
        [Cpv,~] = cpsd(detrend(pp),detrend(vv),hanning(nfft(2)),round(nfft(2)*0.25,0),nfft(2),fs);
        
        %Compute wave energy flux
        zp = AQD.Gatts.transducer_offset_from_bottom; %height of pressure sensor [m]
        df = F(3)-F(2);
        omega = 2*pi.*F;
        k = qkhf(omega,h)./h;
        kh = k*h;
        coshkhz = cosh(k*zp);
        coshkh = cosh(k*h);
        attn = coshkhz./coshkh;                     %transfer function for surface elevation scaling
        attn(attn<0.2) = 0.2;                       %limiter for transfer function
        Spp = Cpp./(attn.^2);                       %surface elevation spectrum
        amp = sqrt(2*Spp.*df);                      %Wave amp = sqrt(2*S(fi)*df)
        E = 0.5*rho*g*amp.^2;                       %Wave energy [Lowe et al., 2005]
        Cg = 0.5*(1+((2*kh)./(sinh(2*kh)))).*(omega./k); %Wave group velocity [Lowe et al., 2005]
        F_aqd(i) = trapz(E(2:end).*Cg(2:end));      %Total Wave energy flux (W/m)
        
        %Compute wave direction (theta)
        Cpu = Cpu.*conj(Cpu);
        Cpv = Cpv.*conj(Cpv);
        Dir=57.296*atan2(Cpu,Cpv);
        Dir=mod(Dir+180,360)-180;                   %Wave direction; 0 deg is north
        DirPeak = Dir(Spp == max(Spp));             %Peak wave direction
        theta_aqd(i) = DirPeak;
    end
end

%Compute wave energy dissipation rate
deltaX = 58.139; %separation distance between C15 and HR ADCP in m
D = abs((F_rbr-F_aqd)./(deltaX*cosd(theta_aqd))); %Total wave energy dissipation rate (W/m^2) [Huang et al., 2012]
%%%CHECK THIS, should this be eq 19 in Lowe??

%Plot figure
ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[800 350   750   400]);
cc = brewermap(2,'Set1');
sp(1) = subplot(211);
pl(1) = plot(time,F_rbr,'linewidth',1.5,'color',cc(1,:));hold on;
pl(2) = plot(time,F_aqd,'linewidth',1.5,'color',cc(2,:));
leg = legend(pl,instName);
set(leg,'position',[0.79 0.825 0.05 0.05])

sp(2) = subplot(212);
plot(time,D,'k','linewidth',1.5)

%Labeling
set(sp,'xlim',[min(time) max(time)])
set(sp(1),'xticklabel',[])
datetick(sp(2),'x','HH:MM','keepticks','keeplimits')
xlabel(sp(2),['\bf\itTime on ' datestr(time(1),'dd-mmm-yyyy')])
ylabel(sp(1),'$F \quad (W/m)$','interpreter','latex')
ylabel(sp(2),'$D \quad (W/m^2)$','interpreter','latex')

%Positioning
set(sp(1),'position',[0.1 0.56 0.8 0.35],...
    'ylim',[50 350],'ytick',50:100:350)
set(sp(2),'position',[0.1 0.13 0.8 0.35])

%Save fig
prettyfigures('text',11,'labels',12,'box',1,'tlength',[0.005 0.005])
set(ff1,'units','inches');
pos = get(ff1,'Position');
set(ff1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(ff1,[figPath 'MKK_C15_HRaqd_WaveEnergyFlux'],'-dpdf','-r0')
