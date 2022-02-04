%Test to check data output from turbulence v mesh size models
clear, close all
path = 'c:\Users\user\Documents\Models\Paper1_2DwaveFlume\ModelRuns\HR_Domain\Scenarios\TEST2\MKKmodel\postProcessing\bathySample\surface\0.25\';
filename = 'epsilon_sampledSurface_bathymetrySample.raw';
fid = fopen([path filename]);
data = textscan(fid,'%n%n%n%n','headerlines',2);
x = data{1};y = data{2};z = data{3};epsilon = data{4};
%There are some negative epsilon values for some reason
epid = find(epsilon>0);
epsilon = epsilon(epid);
x = x(epid);
y = y(epid);
z = z(epid);

[B,I] = sort(z);
figure
plot(x(I),z(I),'.')

figure
plot3(x(I),y(I),z(I),'o')

bins = -5:0.1:-0.5;
%Bin by z-axis
[~,~,loc] = histcounts(z,bins);
lz = find(loc>0);
epsilon_m = accumarray(loc(lz),epsilon(lz))./accumarray(loc(lz),1);
z_m = accumarray(loc(lz),z(lz))./accumarray(loc(lz),1);
figure
plot(epsilon_m,z_m)