% Model canopy density (a) using the Henderson et al. (2017) model
% L = Tw/(2*pi*Tf), where Tf = 2/(Cd*a*ubar) assuming values of the drag
% coefficient Cd and literature values for infragravity waves across coral
% reefs. 
%
% This script is part of the planning phase of my second Mendenhall
% manuscript to determine what density of coral restorations we should be
% aiming for in our numerical models. 
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
figDir = 'c:\Users\bknorris\Documents\Writing\Papers\Paper2_OptimizingRestoration\Planning\';
Cd = [0.5 1.5]; %Osorio-Cano (2019)
Tw = [60 120 240 360]; %infragravity periods (25 - 250 s)
Ub = linspace(0.05,0.5,100); %Wave velocities 

%% Calculate maximum canopy density (a) +/- 10% 
%If we assume L0 = 1, dissipation (X) should be maximized
%Hence we can define L0 = 0.9 for 90% and L0 = 1.1 for 110% of maximum
%value
a_90 = zeros(length(Cd),length(Tw),length(Ub));
a_110 = zeros(length(Cd),length(Tw),length(Ub));
a_max = zeros(length(Cd),length(Tw),length(Ub));
for i = 1:length(Cd)
    for j = 1:length(Ub)
        for k = 1:length(Tw)
            a_max(i,j,k) = (4*pi)/(Cd(i)*Ub(j)*Tw(k)); %Rearranged Eq 18 from Henderson et al. (2017) if Lambda_0 = 1
            a_90(i,j,k) = (0.9*4*pi)/(Cd(i)*Ub(j)*Tw(k));
            a_110(i,j,k) = (1.1*4*pi)/(Cd(i)*Ub(j)*Tw(k));
        end
    end
end
ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[100 80   800   400]);

%FOR CD = 1
sp(1) = subplot(121);
cc = brewermap(length(Tw),'Set1');
p = zeros(length(Tw),1);
for i = 1:length(Tw)
    patch([Ub fliplr(Ub)],...
        [squeeze(a_90(1,1:length(Ub),i)) fliplr(squeeze(a_110(1,1:length(Ub),i)))],...
        cc(i,:),'facealpha',0.25,'edgecolor','none');hold on
    p(i) = plot(Ub,squeeze(a_max(1,1:length(Ub),i)),...
        '-','color',cc(i,:),...
        'markerfacecolor','w',...
        'linewidth',1.2);
end

%FOR CD = 2
sp(2) = subplot(122);
cc = brewermap(length(Tw),'Set1');
p = zeros(length(Tw),1);
legtext = cell(length(Tw),1);
for i = 1:length(Tw)
    patch([Ub fliplr(Ub)],...
        [squeeze(a_90(2,1:length(Ub),i)) fliplr(squeeze(a_110(2,1:length(Ub),i)))],...
        cc(i,:),'facealpha',0.25,'edgecolor','none');hold on
    p(i) = plot(Ub,squeeze(a_max(2,1:length(Ub),i)),...
        '-','color',cc(i,:),...
        'markerfacecolor','w',...
        'linewidth',1.2);
    legtext{i} = ['$T_w = ' sprintf('%.0f',Tw(i)) '\ \mathrm{s}$'];
end
set(sp,'xlim',[0 0.55],'ylim',[0 10])
set(sp(1),'position',[0.1 0.12 0.4 0.8])
set(sp(2),'yticklabel',[],'position',[0.55 0.12 0.4 0.8])
xlabel(sp(1),'$u_b\ \mathrm{(m/s)}$','interpreter','latex')
xlabel(sp(2),'$u_b\ \mathrm{(m/s)}$','interpreter','latex')
ylabel(sp(1),'$a\ \mathrm{(m^{-1})}$','interpreter','latex')
title(sp(1),'$C_D = 0.5$','interpreter','latex')
title(sp(2),'$C_D = 1.5$','interpreter','latex')
prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.008 0.008])
legend(p,legtext,'interpreter','latex','box','off')
% print(ff,[figDir 'Henderson_a_withRanges'],'-dpng','-r0')

%% Other code snippets that might be useful later 

% % Define parameters (comments are literature references)
% Cd = [1 2]; %Osorio-Cano (2019)
% Tw = [25 50 75 100 125 150 200 225 250]; %infragravity periods (25 - 250 s)
% Ub = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5]; %Wave velocities 
% %(estimated from Sous et al. 2020; Pomeroy et al. 2018; Van Dongeren et al. 2013)

% ff1 = figure(1);
% cc = brewermap(length(Tw),'Set1');
% p = zeros(length(Tw),1);
% legtext = cell(length(Tw),1);
% for i = 1:9
%     p(i) = plot(Ub,squeeze(a(1,1:length(Ub),i)),...
%         '-o','color',cc(i,:),...
%         'markerfacecolor','w',...
%         'linewidth',1.2);hold on
%     legtext{i} = ['$T_w = ' sprintf('%.0f',Tw(i)) '\ \mathrm{s}$'];
% end
% set(gca,'xlim',[0 0.55])
% xlabel('$u_b\ \mathrm{(m/s)}$','interpreter','latex')
% ylabel('$a\ \mathrm{(m^{-1})}$','interpreter','latex')
% title('Estimated Canopy Density with C_D = 1')
% prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.008 0.008])
% legend(p,legtext,'interpreter','latex','box','off')
% 
% set(ff1,'units','inches','renderer','painters');
% pos = get(ff1,'Position');
% set(ff1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff1,[figDir 'Henderson_max_a_Cd_1'],'-dpdf','-r0')
% 
% ff2 = figure(2);
% cc = brewermap(length(Tw),'Set1');
% p = zeros(length(Tw),1);
% legtext = cell(length(Tw),1);
% for i = 1:9
%     p(i) = plot(Ub,squeeze(a(2,1:length(Ub),i)),...
%         '-o','color',cc(i,:),...
%         'markerfacecolor','w',...
%         'linewidth',1.2);hold on
%     legtext{i} = ['$T_w = ' sprintf('%.0f',Tw(i)) '\ \mathrm{s}$'];
% end
% set(gca,'xlim',[0 0.55])
% xlabel('$u_b\ \mathrm{(m/s)}$','interpreter','latex')
% ylabel('$a\ \mathrm{(m^{-1})}$','interpreter','latex')
% title('Estimated Canopy Density with C_D = 2')
% prettyfigures('text',11,'labels',12,'box',1,'tickdir','out','tlength',[0.008 0.008])
% legend(p,legtext,'interpreter','latex','box','off')
% 
% set(ff2,'units','inches','renderer','painters');
% pos = get(ff2,'Position');
% set(ff2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff2,[figDir 'Henderson_max_a_Cd_2'],'-dpdf','-r0')

% %% Calculate Henderson model based on parameters
% c = 1;
% a = squeeze(a_max(1,:,1));
% comb = length(Cd)*length(a)*length(Tw)*length(Ub); %combinations of parameters
% L0 = zeros(comb,1);
% X = zeros(comb,1);
% for i = 1:length(Cd)
%     for j = 1:length(Ub)
%         for k = 1:length(Tw)
%             for l = 1:length(a)
%                 Lambda0 = (Cd(i)*a(l)*Ub(j)*Tw(k))/(4*pi);
%                 L0(c) = Lambda0;
%                 X(c) = (((1+4*(Lambda0^2))^(1/2)-1)^(3/2))/((2^(2/3))*(Lambda0^2));
%                 c = c+1;
%             end
%         end
%     end
% end
% figure
% loglog(L0,X,'.')