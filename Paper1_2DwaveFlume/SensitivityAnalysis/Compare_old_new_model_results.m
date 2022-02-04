% Plot free surface spectra, wave energy flux, and total mean TKE
% dissipation to run the sensitivity analysis on model data. Ideally, the
% "new" and "old" models should be nearly identical. 
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

%Define data paths
workingDir = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\DataAnalysis\SensitivityAnalysis\';

modelInfo = 'HR_Domain_model_toRun.csv';

%Load modelInfo from the CSV file and create a structure from the data
fid = fopen([workingDir modelInfo]);
header = textscan(fgetl(fid),repmat('%s',1,7),'delimiter',',');
data = textscan(fid,repmat('%f',1,7),'delimiter',',');
modelInfo = cell2struct(data,[header{:}],2);

%Load "old" models (models that used StokesI wave generation)
oldModels = dir([workingDir '*_old*']);
oldModels = oldModels(~ismember({oldModels.name},{'.','..'}));
oldModels = {oldModels.name};

%Load "new" models (models that used StokesV and cnoidal wave generation)
newModels = dir([workingDir '*_new*']);
newModels = newModels(~ismember({newModels.name},{'.','..'}));
newModels = {newModels.name};

%Loop through model scenarios and plot comparisons
for i = 1:length(oldModels)
    if contains(oldModels{i},'freeSurf')
        load([workingDir oldModels{i}]);
        old = data;clear data
        load([workingDir newModels{i}]);
        new = data;clear data
        modelName = regexprep(oldModels{i},'(_freeSurf_).*','');
        modelName = regexprep(modelName,'_',' ');
        
        ff = figure(1);
        set(ff,'PaperOrientation','landscape',...
            'position',[600 450   1250   400]);

        %Calculate wave spectra @ WG2 (start of hi-res patches), see 
        % c:\Users\user\Documents\Models\Figures\Paper1_2DwaveFlume\HR_LR_calibration_domains_with_annotations_V3.pdf
        [oldEta,F] = pwelch(old.WG2.eta,[],[],[],8);
        [newEta,~] = pwelch(new.WG2.eta,[],[],[],8);
        sp(1) = subplot(131);
        pp1 = plot(F,oldEta,'-k','linewidth',1.5);hold on
        pp2 = plot(F,newEta,'-r','linewidth',1.5);hold on
        set(gca,'yscale','log',...
            'xlim',[0 1])
        leg = legend([pp1 pp2],{'StokesI';'StokesV/cnoidal'});
        xlabel('$f \quad (Hz)$','interpreter','latex')
        ylabel('$S_{\eta} \quad (m^2/Hz)$','interpreter','latex')
        title('Free surface spectra')
        
        %Calculate wave energy flux and dissipation
        scenarioNumber = regexp(modelName,'\ (.*)','tokens');
        modelID = find(strcmp(string(modelInfo.scenarioNumber),scenarioNumber{:}));
        h = modelInfo.h(modelID);
        
        %"Old" models
        [WG1,~] = pwelch(old.WG1.eta,[],[],[],8);
        [WG2,F] = pwelch(old.WG2.eta,[],[],[],8);
        
        k=wavek(F,h);
        kh = k*h;
        c = sqrt(9.81.*tanh(kh))./sqrt(k);
        n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
        Cg = n.*c;
        F1 = trapz(1025*9.81*WG1(2:end).*Cg(2:end));
        F2 = trapz(1025*9.81*WG2(2:end).*Cg(2:end));
        Dold = (F1-F2)/(63-31); %wave energy dissipation between WG1 and WG2 (Huang et al. 2012)
        WEFold = 1025*9.81*WG2(2:end).*Cg(2:end);
        
        %"New" models
        [WG1,~] = pwelch(new.WG1.eta,[],[],[],8);
        [WG2,F] = pwelch(new.WG2.eta,[],[],[],8);
        
        k=wavek(F,h);
        kh = k*h;
        c = sqrt(9.81.*tanh(kh))./sqrt(k);
        n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
        Cg = n.*c;
        F1 = trapz(1025*9.81*WG1(2:end).*Cg(2:end));
        F2 = trapz(1025*9.81*WG2(2:end).*Cg(2:end));
        Dnew = (F1-F2)/(63-31); %wave energy dissipation between WG1 and WG2 (Huang et al. 2012)
        WEFnew = 1025*9.81*WG2(2:end).*Cg(2:end);

        sp(2) = subplot(132);
        pp3 = plot(F(2:end),WEFold./1000,'-k','linewidth',1.5);hold on
        pp4 = plot(F(2:end),WEFnew./1000,'-r','linewidth',1.5);
        set(gca,'yscale','log',...
            'xlim',[0 1])
        leg = legend([pp1 pp2],{'StokesI';'StokesV/cnoidal'});
        xlabel('$f \quad (Hz)$','interpreter','latex')
        ylabel('$F_0 \quad (KW/m)$','interpreter','latex')
        title('Wave energy flux')
        
        sp(3) = subplot(133);
        pp5 = bar(categorical({'StokesI','StokesV/cnoidal'}),[Dold Dnew]);
        ylabel('$D \quad (m^2/s^3)$','interpreter','latex')
        title('Total Wave Dissipation')
        
        suptitle(modelName)
        
        RMSE = sqrt(mean((Dold-Dnew).^2));
        pctErr = ((Dold-Dnew)/Dnew)*100;
        fprintf('RMSE between stokesI and stokesV total wave dissipation = %0.2f\n',RMSE)
        fprintf('Percent diff. between stokesI and stokesV total wave dissipation = %0.2f\n\n',pctErr)
    elseif contains(oldModels{i},'rawData')
        load([workingDir oldModels{i}]);
        old = epsilon;clear epsilon
        load([workingDir newModels{i}]);
        new = epsilon;clear epsilon
        modelName = regexprep(oldModels{i},'(rawData).*','');
        modelName = regexprep(modelName,'_',' ');
        
        ff = figure(2);
        set(ff,'PaperOrientation','landscape',...
            'position',[600 450   1250   400]);
        %Plot turbulence data
        %Normalize turbulence by incoming wave energy dissipation

        
        %"Old" models
        epsAvg = mean(old.epsilon,2);
        x = old.x(:,1)-30;
        z = old.z(:,1);
        epsOld = epsAvg./Dold;
        
        nBins = 100;
        xBins = linspace(min(x),max(x),nBins);
        zBins = linspace(min(z),max(z),nBins);
        
        ix = discretize(x, xBins);
        iz = discretize(z, zBins);
        idx = sub2ind([nBins nBins], ix, iz);
        epsg = accumarray(idx, epsOld, [nBins*nBins 1], @(x) mean(x), NaN);
        epsg = reshape(epsg, nBins, nBins);
        
        sp(1) = subplot(131);
        p = pcolor(xBins,zBins,epsg');
        shading interp
        hold on
        set(gca,'ydir','normal')
        cc = brewermap(nBins,'*RdYlBu');
        colormap(cc)
        caxis([1e-9 1e-5])
        cb1 = colorbar(sp(1),'location','southoutside');
        
        %"New" models
        epsAvg = mean(new.epsilon,2);
        x = new.x(:,1)-30;
        z = new.z(:,1);
        epsNew = epsAvg./Dnew;
        
        xBins = linspace(min(x),max(x),nBins);
        zBins = linspace(min(z),max(z),nBins);
        
        ix = discretize(x, xBins);
        iz = discretize(z, zBins);
        idx = sub2ind([nBins nBins], ix, iz);
        epsg = accumarray(idx, epsNew, [nBins*nBins 1], @(x) mean(x), NaN);
        epsg = reshape(epsg, nBins, nBins);
        
        sp(2) = subplot(132);
        p = pcolor(xBins,zBins,epsg');
        shading interp
        hold on
        set(gca,'ydir','normal')
        cc = brewermap(nBins,'*RdYlBu');
        colormap(cc)
        caxis([1e-9 1e-5])
        cb2 = colorbar(sp(2),'location','southoutside');
        
        %Differenced models
        sp(3) = subplot(133);
        diffEps = epsOld-epsNew;
        xBins = linspace(min(x),max(x),nBins);
        zBins = linspace(min(z),max(z),nBins);
        
        ix = discretize(x, xBins);
        iz = discretize(z, zBins);
        idx = sub2ind([nBins nBins], ix, iz);
        epsg = accumarray(idx, diffEps, [nBins*nBins 1], @(x) mean(x), NaN);
        epsg = reshape(epsg, nBins, nBins);
        
        sp(3) = subplot(133);
        p = pcolor(xBins,zBins,epsg');
        shading interp
        hold on
        set(gca,'ydir','normal')
        cc = brewermap(nBins,'*RdYlBu');
        colormap(cc)
        caxis([1e-9 1e-5])
        cb3 = colorbar(sp(3),'location','southoutside');
        
        title(sp(1),'StokesI Model')
        title(sp(2),'StokesV/cnoidal Model')
        title(sp(3),'Difference')
        ylabel(sp(1),'z (m)')
        ylabel(sp(2),'z (m)')
        ylabel(sp(3),'z (m)')
        xlabel(sp(1),'Across-shore Distance (m)')
        xlabel(sp(1),'Across-shore Distance (m)')
        xlabel(sp(3),'Across-shore Distance (m)')
        xlabel(cb1,'$\overline{\epsilon}/\overline{\epsilon_0} \quad (-)$','interpreter','latex')
        xlabel(cb2,'$\overline{\epsilon}/\overline{\epsilon_0} \quad (-)$','interpreter','latex')
        xlabel(cb3,'$\overline{\epsilon}/\overline{\epsilon_0} \quad (-)$','interpreter','latex')
        suptitle(modelName)
        
        epsTotOld = trapz(epsOld);
        epsTotNew = trapz(epsNew);
        RMSE = sqrt(mean((epsTotOld-epsTotNew).^2));
        pctErr = ((epsTotOld-epsTotNew)/epsTotNew)*100;
        fprintf('RMSE between stokesI and stokesV total TKE dissipation = %0.3f\n',RMSE)
        fprintf('Percent diff. between stokesI and stokesV total TKE dissipation = %0.2f\n\n',pctErr)
    end
    prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.01 0.01])
    pause()
    close all
end



