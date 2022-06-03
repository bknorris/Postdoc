
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% RUNUP BASED ON RBF INTERPOLATION FROM X-BEACH SIMULATIONS ON CORAL REEFS  %%%%%%%%%%%%%%%%%%%% 
%%%%%% RUNUP = RBF{'WL'; 'Hm0'; 'Tp'; 'Rslope';'Bslope';'ReefW'; 'Cf'} FOR SEMI-INFINITE SLOPE %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
% SANTANDER, MAY 2017
% CODED BY SURF & SURGE RESEARCH GROUP, UNIVERSIDAD DE CANTABRIA
% XBEACH SIMULATIONS BY DELTARES AND USGS SANTA CRUZ (PEARSON 2016)

clear all; close all; clc;

addpath([pwd  filesep 'RBF' filesep])

%% LOAD 
% DATA MATRIX
% Structure: ['WL', 'Hm0', 'Tp', 'Reef_slope','Beach_slope','ReefW','Cf']

load('Input_data.mat'); % waves and reef conditions
data_v=Input_data;
%data_v = [2.73,2.1,15,0.5,0.2,100,0.04;-0.73,2.1,15,0.5,0.2,100,0.04]%;0.61,2.1,15,0.5,0.2,100,0.04;0.98,1.6,13,0.5,0.2,100,0.04];

data1 = data_v(:,[1,4:7]);

%nans removed
posn = find(isnan(data1(:,end))==1);
data1(posn,:) = [];

[N,dim]=size(data1);

%min and max
load([pwd filesep 'RBF_coefficients' filesep 'Min_from_simulations.mat'])
load([pwd filesep 'RBF_coefficients' filesep 'Max_from_simulations.mat'])
minimum1 = minimum([1,4:7]);
maximum1 = maximum([1,4:7]);

%% INPUT DATA QUALITY CHECK
hs_lo2=data_v(:,2)./(9.81/(2*pi).*data_v(:,3).^2);

out=zeros(size(data_v,1),1);

for i=1:size(data_v,1)
    for j=1:size(data_v,2)
        if data_v(i,j)<minimum(j) || data_v(i,j)>maximum(j)
            out(i)=max(out(i),1);
        else
            if hs_lo2(i)<0.005-0.0001 || hs_lo2(i)>0.05+0.0001
                out(i)=max(out(i),1);
            else
                out(i)=max(out(i),0);
            end
        end
    end
end

s=find(out==1); %rows out of bound


%% INTERPOLATION RESULTS

if isempty(s)
    % data normalization
    datos_n2=zeros(N,dim);
    for i=1:size(data1,2)
        datos_n2(:,i)=(data1(:,i)-minimum1(i))/(maximum1(i)-minimum1(i));
    end
    
    ncases=15; %Number of wave conditions simulated
    runup=nan(size(data_v,1),ncases);
    
    %interpolations
    %rbf
    for k=1:size(data_v,1)
        for i=1:ncases
            
            load([pwd filesep 'RBF_coefficients' filesep 'Coeffs_Runup_Xbeach_test' num2str(i) '.mat'])
            runup(k,i) = rbfinterp_modificado(datos_n2(k,1:end)',coeffsave);
            
        end
    end
    
    
    hs=[1 1 1 2 2 2 3 3 3 4 4 4 5 5 5]; %Wave heigth conditions used in RBF
    hs_lo=[0.005 0.025 0.05 0.005 0.025 0.05 0.005 0.025 0.05 0.005 0.025 0.05 0.005 0.025 0.05 ]; %Wave conditions used in RBF
    
    %linear
    for j=1:size(data_v,1)
        
        x=hs; y=hs_lo; z=runup(j,:);
        hs_e=data_v(j,2);
        tp_e=data_v(j,3);
        lo=9.81/(2*pi).*tp_e.^2;
        hs_lo_e=hs_e/lo;
        x_e=hs_e; y_e=hs_lo_e;
        vq = griddata(x,y,z,x_e,y_e);
        
        RU(j)=vq;
        
    end
    
    %  figure;
    %  mesh(x_e,y_e,vq), hold on, plot3(x,y,z,'o'), hold off
    %
    RMSE=RU*0.1653; %0.1653 = mean SI from kfold valdation
    
    disp('  RUNUP(m)      RMSE(m) ')
    disp([RU',RMSE'])
    
    save('Output.mat','RU','RMSE')
    
else
    disp(['ERROR: The following rows of data are out of bounds: ' num2str(s')])
end
