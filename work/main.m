% Wavelet-approximation by LSE
%
% Version: one-level. Date: 31.01.2025
close all; clc;

% для 400
% -11 -8    20 -- оптимально по rmse
% -9  -9    100 -- интересно по уровням
% -12 -7    20 -- наибольший размах

% для 1000
% -9  -9  100 -- оптимально по rmse
% -10 -8  40 -- несколько уровней


% Wavelet settings
j_min = -20;
j_max = -15;       % resolution level 
k_max = 40;      % nr. of wavelets on t-axis



Filename_ggm = '../data/XGM2019_400.dat';
fid        =  fopen(Filename_ggm,'r');              
ggm400        =  textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',44);  
fclose(fid);

Filename_ggm = '../data/XGM2019_1000.dat';
fid        =  fopen(Filename_ggm,'r');              
ggm1000        =  textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',44);  
fclose(fid);

TimeGPS = ggm400{1};

DG3 = 10^(-5)*(ggm1000{5}-ggm400{5});
DG1 = -deg2rad(ggm1000{6}-ggm400{6}/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});
DG2 = -deg2rad(ggm1000{7}-ggm400{7}/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});

start   = 200;
fin     = 1478;

% Interval & Input function
TimeArray = TimeGPS(start:fin) - TimeGPS(start);
DG1        = DG1(start:fin); 
DG2        = DG2(start:fin); 
DG3        = DG3(start:fin); 
Time_fin  = TimeArray(end);


% Matrix of wavelet values

Psi = zeros( length(TimeArray), (j_max - j_min + 1) * k_max );
for j = j_min : j_max
    for k = 1 : k_max
       dT  = (2^j)*Time_fin/k_max;  % step of wavelet grid
       t_k = k * dT;                % knot in grid
       Psi(:,(j - j_min) * k_max + k) = wavel_trf(j,t_k,TimeArray);   
    end
end
OneDim = 3;
if OneDim ~= 1
    Psi = [Psi, zeros(size(Psi)), zeros(size(Psi));
           zeros(size(Psi)),  Psi, zeros(size(Psi));
           zeros(size(Psi)), zeros(size(Psi)), Psi;];
    
    DG = [DG1; DG2; DG3];
else
    DG = DG3;
end

% LS estimation
WCoeff = Psi \ DG;
disp(['Nr. of wavelet-coef: ',num2str(max(size(WCoeff)))])

disp(['Analysing SVD...'])
s=svd(Psi);
if s(1)/s(end)>10^16
disp(['Warning: matrix is ill-conditioned'])    
end
disp(['Condition nr: ',num2str(s(1)/s(end))])


% Wavelet-reconstruction
DG_est = Psi * WCoeff;

DG_est1 = DG_est(1:end/3, :);
DG_est2 = DG_est(end/3+1:2*end/3, :);
DG_est3 = DG_est(2*end/3+1:end, :);
% DG_est3 = DG_est;

figure(1)
plot(TimeArray,DG1)
hold on;
plot(TimeArray, DG_est1,'r')
legend('true','estimate')
grid on;
title(['Reconstructed input function. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
xlabel('Time(s)')

figure(2)
plot(TimeArray,DG2)
hold on;
plot(TimeArray, DG_est2,'r')
legend('true','estimate')
grid on;
title(['Reconstructed input function. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
xlabel('Time(s)')

figure(3)
plot(TimeArray,DG3)
hold on;
plot(TimeArray, DG_est3,'r')
legend('true','estimate')
grid on;
title(['Reconstructed input function. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
xlabel('Time(s)')

delta = j_max - j_min +1;

WCoeff1 = WCoeff(1:length(WCoeff)/3);
figure(4)
for i = 1:delta
    subplot(delta,1,i)
    WCoeffSub = WCoeff1((i-1)/delta*end+1:i/delta*end);
    plot(WCoeffSub); hold on;
    plot(WCoeffSub,'.')
    title(['\Delta g_1 coeffiсients, coeffiсientslevel=',num2str(j_min+i-1)])
    xlabel('Number of coefficient')
end

WCoeff2 = WCoeff(length(WCoeff)/3+1:2*length(WCoeff)/3);
figure(5)
for i = 1:delta
    subplot(delta,1,i)
    WCoeffSub = WCoeff2((i-1)/delta*end+1:i/delta*end);
    plot(WCoeffSub); hold on;
    plot(WCoeffSub,'.')
    title(['\Delta g_2 coeffiсients, level=',num2str(j_min+i-1)])
    xlabel('Number of coefficient')
end

WCoeff3 = WCoeff(2*length(WCoeff)/3+1:end);
figure(6)
for i = 1:delta
    subplot(delta,1,i)
    WCoeffSub = WCoeff3((i-1)/delta*end+1:i/delta*end);
    plot(WCoeffSub); hold on;
    plot(WCoeffSub,'.')
    title(['\Delta g_3 coeffiсients, level=',num2str(j_min+i-1)])
    xlabel('Number of coefficient')
end
