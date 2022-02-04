%%%COMPARE INSTRUMENT FIELDS TEST TEST TEST%%%
clear
load('MKK_Velocities_test.mat')
load('MKK_HR-Aquadopp_zuv.mat')
load('MKK_model_depthM.mat')

%TEST CODE FOR INSTRUMENT DATA
[m,n]=size(dat.U);
Uwrms_i = zeros(m,1);
Wwrms_i = zeros(m,1);
for i = 1:m
    U = dat.U(i,:);
    V = dat.V(i,:);
    W = dat.W(i,:);
    nwin = length(U);
    
    %Bridge short gaps in time series if NaNs are present
    if any([isnan(U) isnan(V) isnan(W)])
        disp('Found NaNs in Velocity... interpolating with cmgbridge')
        fprintf('Bin : %0.f\n',i)
        U = cmgbridge(U,100,100,1000);
        V = cmgbridge(V,100,100,1000);
        W = cmgbridge(W,100,100,1000);
    end
    
    %Calculate RMS velocities
    Ec = (1/nwin)*sum(U);Nc = (1/nwin)*sum(V);Wc = (1/nwin)*sum(W);
    Ewrms = sqrt((1/nwin)*sum((U-Ec).^2));
    Nwrms = sqrt((1/nwin)*sum((V-Nc).^2));
    Wwrms_i(i,:) = sqrt((1/nwin)*sum((W-Wc).^2));
    Uwrms_i(i,:) = sqrt(Ewrms^2+Nwrms^2);
end
figure(1)
plot(Uwrms_i,zuv,'+b'), hold on
plot(Wwrms_i,zuv,'+r')

clear m n U V W
[m,n]=size(dat.Um);
Uwrms_m = zeros(m,1);
Wwrms_m = zeros(m,1);
for i = 1:m
    U = dat.Um(i,:);U = U(~isnan(U));
    V = dat.Vm(i,:);V = V(~isnan(V));
    W = dat.Wm(i,:);W = W(~isnan(W));
    nwin = length(U);
    
    %Bridge short gaps in time series if NaNs are present
    if any([isnan(U) isnan(V) isnan(W)])
        disp('Found NaNs in Velocity... interpolating with cmgbridge')
        fprintf('Bin : %0.f\n',i)
        U = cmgbridge(U,100,100,1000);
        V = cmgbridge(V,100,100,1000);
        W = cmgbridge(W,100,100,1000);
    end
    
    %Calculate RMS velocities
    Ec = (1/nwin)*sum(U);Nc = (1/nwin)*sum(V);Wc = (1/nwin)*sum(W);
    Ewrms = sqrt((1/nwin)*sum((U-Ec).^2));
    Nwrms = sqrt((1/nwin)*sum((V-Nc).^2));
    Wwrms_m(i,:) = sqrt((1/nwin)*sum((W-Wc).^2));
    Uwrms_m(i,:) = sqrt(Ewrms^2+Nwrms^2);
end

plot(Uwrms_m,depthM,'-b')
plot(Wwrms_m,depthM,'-r')