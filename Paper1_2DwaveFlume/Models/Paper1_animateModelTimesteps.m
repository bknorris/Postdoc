% Utility to plot all timesteps of a particular model
%
%
% Updates: added a loop to plot all processed models to look for oddities
% in the data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all
modelPath = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\DataAnalysis\HR_Domain\';
figPath = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\Figures\ModelAnimations\';

modelNames = dir([modelPath '*rawData_V2*']);
modelNames = modelNames(~ismember({modelNames.name},{'.','..'}));
modelNames = {modelNames.name};
freeSurf = dir([modelPath '*freeSurf_V2*']);
freeSurf = freeSurf(~ismember({freeSurf.name},{'.','..'}));
freeSurf = {freeSurf.name};
for i = 1:length(modelNames)
    fprintf('Loading file %0.0f of %0.0f: %s\n',i,length(modelNames),modelNames{i})
    model = regexp(modelNames{i},'(?<=_)(.*?)(?=_)','match');
    if contains(modelPath,'HR')
        videoTitle = ['HR Scenario ' model{1}];
        videoFile = ['HR_Domain_Scenario_' model{1}];
    elseif contains(modelPath,'LR')
        videoTitle = ['LR Scenario ' model{1}];
        videoFile = ['LR_Domain_Scenario_' model{1}];
    end
    if isfile([figPath videoFile '.avi'])
        disp('The animation has already been created!')
        continue
    else
        data = load([modelPath modelNames{i}],'time','U','epsilon');
        FS = load([modelPath freeSurf{i}]);
    end
    textX = 0.8;
    textY = 0.1;
    
    ff = figure(1);
    set(ff,'PaperOrientation','landscape',...
        'position',[500 80   650   800]);
    video = VideoWriter([figPath videoFile '.avi']);
    video.FrameRate = 10;
    open(video)
    timeSteps = min([max(size(data.time)) min(size(data.U.x)) min(size(data.epsilon.x))]);
    for j = 1:timeSteps
        title(videoTitle)
        sp(1) = subplot(311);
        x = FS.waveGauges(2:end)-30;
        waveGauges = fieldnames(FS);
        eta = zeros(length(waveGauges(1:end-2)),1);
        for k = 1:length(waveGauges(1:end-2))
            eta(k) = FS.(waveGauges{k}).eta(j);
        end
        plot(x,eta(2:end),'-o','color',[0 57 255]./255,'linewidth',1.5,'markerfacecolor',[0 57 255]./255,'markersize',5)
        ylabel('$\eta \ \mathrm{(m)}$','interpreter','latex')
        set(sp(1),'ylim',[-2 2],...
            'xticklabel',[],...
            'position',[0.15 0.86 0.7 0.1])
        
        sp(2) = subplot(312);
        x = data.U.x(:,j)-30;
        z = data.U.z(:,j);
        Umag = sqrt((data.U.Ux(:,j).^2)+(data.U.Uy(:,j).^2)+(data.U.Uz(:,j).^2));
        
        nBins = 100;
        xBins = linspace(min(x),max(x),nBins);
        zBins = linspace(min(z),max(z),nBins);
        
        ix = discretize(x, xBins);
        iz = discretize(z, zBins);
        idx = sub2ind([nBins nBins], ix, iz);
        Umagg = accumarray(idx, Umag, [nBins*nBins 1], @(x) mean(x), NaN);
        Umagg = reshape(Umagg, nBins, nBins);
        p1 = pcolor(xBins,zBins,Umagg');
        text(textX,textY,sprintf('t = %0.3f s',data.time(j)),'units','normalized')
        shading interp
        cc = brewermap(nBins,'*RdYlBu');
        colormap(cc)
        caxis([0 1.5])
        cb = colorbar;
        set(cb,'location','eastoutside')
        ylabel(cb,'$\mathrm{Umag \ (m/s)}$','interpreter','latex')
        ylabel('$z \ \mathrm{(m)}$','interpreter','latex')
        set(sp(2),'xticklabel',[],...
            'ydir','normal',...
            'position',[0.15 0.48 0.7 0.35])
        
        sp(3) = subplot(313);
        x = data.epsilon.x(:,j)-30;
        z = data.epsilon.z(:,j);
        eps = data.epsilon.epsilon(:,j);
        nBins = 100;
        xBins = linspace(min(x),max(x),nBins);
        zBins = linspace(min(z),max(z),nBins);
        
        ix = discretize(x, xBins);
        iz = discretize(z, zBins);
        idx = sub2ind([nBins nBins], ix, iz);
        epsg = accumarray(idx, eps, [nBins*nBins 1], @(x) mean(x), NaN);
        epsg = reshape(epsg, nBins, nBins);
        p2 = pcolor(xBins,zBins,epsg');
        text(textX,textY,sprintf('t = %0.3f s',data.time(j)),'units','normalized')
        shading interp
        cc = brewermap(nBins,'*RdYlBu');
        colormap(cc)
        caxis([10E-5 10E-2])
        cb = colorbar;
        set(cb,'location','eastoutside')
        ylabel(cb,'$\varepsilon \ \mathrm{(m^2/s^3)}$','interpreter','latex')
        xlabel('Across shore distance (m)','interpreter','latex')
        ylabel('$z \ \mathrm{(m)}$','interpreter','latex')
        set(sp(3),'ydir','normal','position',[0.15 0.1 0.7 0.35])
        prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.01 0.01])
        set(sp,'xlim',[round(min(x),0) round(max(x),0)],...
            'xtick',round(min(x),0):1:round(max(x),0))
        pause(0.01)
        frame = getframe(gcf);
        writeVideo(video, frame);
    end
    close(video)
    close(ff)
end