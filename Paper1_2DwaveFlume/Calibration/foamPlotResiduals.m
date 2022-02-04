% Plot residuals from an OpenFOAM run
%
% Remember, run the command "foamLog log.interFoam" in Ubuntu before trying 
% this script to create the 'logs' folder in the model directory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close
fdir = 'c:\Users\user\Documents\Models\MKK_Combined\2D\Burst151-158\MKK_CombinedHR_V11\MKKmodel\logs\';
fields = {'p_rgh_2';'p_rghFinalRes_2';'k_0';'kFinalRes_0';'epsilon_0';'epsilonFinalRes_0'}; 
figPath = 'c:\Users\user\Documents\Models\Figures\';
modelName = 'MKK_CombinedHR_V15_burst151-158';

ff = figure(1);
set(ff,'PaperOrientation','landscape',...
    'position',[200 400   1200   450]);
hold on
cc = brewermap(length(fields),'Paired');
p = zeros(size(fields));

for i = 1:length(fields)
    fid = fopen([fdir fields{i}]);
    text = textscan(fid,'%n%n');
    p(i) = plot(text{1},text{2},'-','linewidth',1.5,'color',cc(i,:));
end
leg = legend(p,fields);
set(leg,'position',[0.9 0.35 0.05 0.2])
set(gca,'xlim',[min(text{1}) max(text{1})],...
    'ylim',[10^-10 10^0],...
    'yscale','log',...
    'position',[0.1 0.15 0.75 0.8])
xlabel('Model Time (s)')
ylabel('Residuals')
prettyfigures('text',12,'labels',13,'tickdir','out','tlength',[0.005 0.005])
hold off
% export_fig(ff,[figPath modelName '_residuals'],'-png')
