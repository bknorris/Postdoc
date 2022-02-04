%Compare model derived quantities (e.g., U, epsilon) with AQDP measurements.
%
% Notes: Unlike the MKK_Compare_Model_Inst_Fields scripts, this script only
% works for lines extracted from the model over the same range as the
% Aquadopp profile (underwater). Do not run with 'old' profile data!
%
% Updates:
% 03/15/21: Revised script to plot both mean flow profiles and scatterplot
% model vs observations. 
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all
dataPath = 'c:\Users\user\Documents\Models\DataAnalysis\Calibration\';
figPath = 'c:\Users\user\Documents\Models\Figures\'; %where to save figs
modelFigName = 'HR_V27_LR_V2';
fileName = dir(dataPath);
fileName = fileName(~ismember({fileName.name},{'.','..'}));
fileName = {fileName.name};
load([dataPath fileName{1}]);cal.HR = data;
load([dataPath fileName{2}]);cal.LR = data;clear data

%% Plot the calibrations
ff1 = figure(1);
set(ff1,'PaperOrientation','landscape',...
    'position',[800 350   1000   600]);
cc = brewermap(2,'Set1');
p1 = zeros(2,1);
p2 = zeros(2,1);
marker = {'+','x'};
textLocX = [0.085 0.085 0.085 0.003 0.003 0.003 3E-6 3E-6 3E-6;...
    0.134 0.134 0.134 0.033 0.033 0.033 2.4E-5 2.4E-5 2.4E-5];
textLocY = [0.137 0.14 0.143 0.03 0.0325 0.035 2.3E-5 2.5E-5 2.7E-5;...
    0.103 0.106 0.109 0.0025 0.005 0.0075 3.0E-6 5.0E-6 7E-6];
for i = 1:2
    fn = fieldnames(cal);
    depth_i = cal.(fn{i}).depth_i;
    depth_m = cal.(fn{i}).depth_m;
    idi = find(depth_i>=0,1,'last');
    idm = find(depth_m>=depth_i(idi) & depth_m<=max(depth_i));
    
    %Bin model depths by ADCP velocity bin heights
    [~,~,loc] = histcounts(depth_m(idm),flipud(depth_i(1:idi)));
    
    Uwrms_i = nanmean(cal.(fn{i}).Uwrms_i(1:idi,:),2);
    Uwrms_m = accumarray(loc(:),flipud(cal.(fn{i}).Uwrms_m(idm)))./accumarray(loc(:),1);
    
    Wwrms_i = nanmean(cal.(fn{i}).Wwrms_i(1:idi,:),2);
    Wwrms_m = accumarray(loc(:),flipud(cal.(fn{i}).Wwrms_m(idm)))./accumarray(loc(:),1);
    
    E_i = nanmean(cal.(fn{i}).E_i(1:idi,:),2);    
    E_m = accumarray(loc(:),flipud(cal.(fn{i}).E_m(idm)))./accumarray(loc(:),1);
    depth_m = accumarray(loc(:),flipud(depth_m(idm)))./accumarray(loc(:),1);
    
    %Plot results
    sp(1) = subplot(231);
    p1(i) = plot(Uwrms_i,depth_i(1:idi),marker{i},...
        'markerfacecolor',cc(i,:),...
        'linewidth',1.5,'color',cc(i,:));hold on
    p2(i) = plot(cal.(fn{i}).Uwrms_m,cal.(fn{i}).depth_m,'linewidth',1.5,'color',cc(i,:));
    
    sp(2) = subplot(232);
    plot(Wwrms_i,depth_i(1:idi),marker{i},...
        'markerfacecolor',cc(i,:),...
        'linewidth',1.5,'color',cc(i,:));hold on
    plot(cal.(fn{i}).Wwrms_m,cal.(fn{i}).depth_m,'linewidth',1.5,'color',cc(i,:));
    
    sp(3) = subplot(233);
    p1(i) = plot(E_i,depth_i(1:idi),marker{i},...
        'markerfacecolor',cc(i,:),...
        'linewidth',1.5,'color',cc(i,:));hold on
    p2(i) = plot(cal.(fn{i}).E_m,cal.(fn{i}).depth_m,'linewidth',1.5,'color',cc(i,:));

    %Calculate and plot model-observation scatter plots
    sp(4) = subplot(234);
    px = 1:10;
    py = 1:10;
    plot(Uwrms_i(px),Uwrms_m(py),'.','markersize',12,'color',cc(i,:)); hold on
    
    %Scatter index (SCI) and Bias [from Van Der Westhuysen, 2010: JGRO]
    SCI = sqrt((1/length(Uwrms_m(py)))*sum((Uwrms_m(py)-Uwrms_i(px)).^2))./((1/length(Uwrms_m(py)))*sum(Uwrms_i(px)));
    Bias = (sum(Uwrms_m(py)-Uwrms_i(px)))./((1/length(Uwrms_m(py)))*sum(Uwrms_i(px)));
    
    %Linear regression and R^2
    pf = polyfit(Uwrms_i(px),Uwrms_m(py),1);
    xx = linspace(0.05,0.15,10);
    pv = polyval(pf,xx);
    R = corrcoef(Uwrms_i(px),Uwrms_m(py));
    Rsq = R(1,2).^2;
    plot(xx,pv,'linewidth',1.5,'color',cc(i,:))
    
    %Add SCI, Bias, and R-squared to figure
    text(textLocX(i,1),textLocY(i,1),sprintf('%0.2f',SCI),'color',cc(i,:))
    text(textLocX(i,2),textLocY(i,2),sprintf('%0.2f',Bias),'color',cc(i,:))
    text(textLocX(i,3),textLocY(i,3),sprintf('%0.2f',Rsq),'color',cc(i,:))
    
    sp(5) = subplot(235);
    if i == 1
        px = 5:11;
        py = 5:11;
    else
        px = 2:10;
        py = 2:10;
    end
%     px = 2:length(Wwrms_i);
%     py = 1:length(Wwrms_m);
    plot(Wwrms_i(px),Wwrms_m(py),'.','markersize',12,'color',cc(i,:)); hold on
    
    %Scatter index (SCI) and Bias [from Van Der Westhuysen, 2010: JGRO]
    SCI = sqrt((1/length(Wwrms_m(py)))*sum((Wwrms_m(py)-Wwrms_i(px)).^2))./((1/length(Wwrms_m(py)))*sum(Wwrms_i(px)));
    Bias = (sum(Wwrms_m(py)-Wwrms_i(px)))./((1/length(Wwrms_m(py)))*sum(Wwrms_i(px)));
    
    %Linear regression and R^2
    pf = polyfit(Wwrms_i(px),Wwrms_m(py),1);
    xx = linspace(0,0.05,10);
    pv = polyval(pf,xx);
    R = corrcoef(Wwrms_i(px),Wwrms_m(py));
    Rsq = R(1,2).^2;
    plot(xx,pv,'linewidth',1.5,'color',cc(i,:))
    
    %Add SCI, Bias, and R-squared to figure
    text(textLocX(i,4),textLocY(i,4),sprintf('%0.2f',SCI),'color',cc(i,:))
    text(textLocX(i,5),textLocY(i,5),sprintf('%0.2f',Bias),'color',cc(i,:))
    text(textLocX(i,6),textLocY(i,6),sprintf('%0.2f',Rsq),'color',cc(i,:))
    
    sp(6) = subplot(236);
    if i == 1
        px = 5:13;
        py = 5:13;
    else
     px = 5:14;
    py = 5:14;   
    end
    plot(E_i(px),E_m(py),'.','markersize',12,'color',cc(i,:)); hold on
    
    %Scatter index (SCI) and Bias [from Van Der Westhuysen, 2010: JGRO]
    SCI = sqrt((1/length(E_m(py)))*sum((E_m(py)-E_i(px)).^2))./((1/length(E_m(py)))*sum(E_i(px)));
    Bias = (sum(E_m(py)-E_i(px)))./((1/length(E_m(py)))*sum(E_i(px)));
    
    %Linear regression and R^2
    pf = polyfit(E_i(px),E_m(py),1);
    xx = linspace(1E-6,1E-4,10);
    pv = polyval(pf,xx);
    R = corrcoef(E_i(px),E_m(py));
    Rsq = R(1,2).^2;
    plot(xx,pv,'linewidth',1.5,'color',cc(i,:))
    
    %Add SCI, Bias, and R-squared to figure
    text(textLocX(i,7),textLocY(i,7),sprintf('%0.2f',SCI),'color',cc(i,:))
    text(textLocX(i,8),textLocY(i,8),sprintf('%0.2f',Bias),'color',cc(i,:))
    text(textLocX(i,9),textLocY(i,9),sprintf('%0.2f',Rsq),'color',cc(i,:))
end
%Labeling
xlabel(sp(1),'$U_{w,RMS} \quad \mathrm{(m/s)}$','interpreter','latex')
xlabel(sp(2),'$W_{w,RMS} \quad \mathrm{(m/s)}$','interpreter','latex')
xlabel(sp(3),'$\varepsilon \quad \mathrm{(m^2/s^2)}$','interpreter','latex')
xlabel(sp(4),'$U_{w,RMS} \ \mathrm{Measured \quad (m/s)}$','interpreter','latex')
xlabel(sp(5),'$W_{w,RMS} \ \mathrm{Measured \quad (m/s)}$','interpreter','latex')
xlabel(sp(6),'$\epsilon \ \mathrm{Measured \quad (m^2/s^2)}$','interpreter','latex')

ylabel(sp(1),'$h \quad \mathrm{(m)}$','interpreter','latex')
ylabel(sp(4),'$U_{w,RMS}\ \mathrm{Modeled \quad (m/s)}$','interpreter','latex')
ylabel(sp(5),'$W_{w,RMS} \ \mathrm{Modeled \quad (m/s)}$','interpreter','latex')
ylabel(sp(6),'$\varepsilon \ \mathrm{Modeled \quad (m^2/s^2)}$','interpreter','latex')

%Axes adjustments
set([sp(1) sp(2) sp(3)],'ylim',[0 0.7],'ytick',0:0.1:0.7)
set(sp(1),'xlim',[0 0.15])
set(sp(2),'xlim',[0 0.10])
set(sp(3),'xscale','log','xlim',[10^-7 10^-3])
set(sp(4),'xlim',[0.08 0.15],'ylim',[0.1 0.15])
set(sp(5),'xlim',[0 0.04],'ylim',[0 0.04])
set(sp(6),'xlim',[10^-6 3e-5],'ylim',[10^-6 3e-5])

%Positioning
set(sp(1),'position',[0.1 0.55 0.25 0.42])
set(sp(2),'position',[0.4 0.55 0.25 0.42])
set(sp(3),'position',[0.7 0.55 0.25 0.42])

set(sp(4),'position',[0.12 0.1 0.21 0.32])
set(sp(5),'position',[0.42 0.1 0.21 0.32])
set(sp(6),'position',[0.72 0.1 0.21 0.32])

prettyfigures('text',11,'labels',10)
set(ff1,'units','inches');
pos = get(ff1,'Position');
set(ff1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% export_fig(ff1,[figPath modelFigName '_AQDP_MODEL-Cals_V2'],'-pdf','-nocrop','-nofontswap')
print(ff1,[figPath modelFigName '_AQDP_MODEL-Cals_V2'],'-dpdf','-r0')
