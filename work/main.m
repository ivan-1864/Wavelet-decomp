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


% для сложного условия
% -8 -7 80
% -8 -8 80

% Wavelet settings
j_min = -8;
j_max = -6;       % resolution level 
k_max = 30;      % nr. of wavelets on t-axis



Filename_ggm = '../data/XGM2019_400.dat';
fid        =  fopen(Filename_ggm,'r');              
ggm400        =  textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',44);  
fclose(fid);

Filename_ggm = '../data/XGM2019_1000.dat';
fid        =  fopen(Filename_ggm,'r');              
ggm1000        =  textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',44);  
fclose(fid);

TimeGPS = ggm400{1};

DG3_res = 10^(-5)*(ggm1000{5}-ggm400{5});
DG1_res = -deg2rad((ggm1000{6}-ggm400{6})/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});
DG2_res = -deg2rad((ggm1000{7}-ggm400{7})/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});

DG3_low = 10^(-5)*(ggm400{5});
DG1_low = -deg2rad(ggm400{6}/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});
DG2_low = -deg2rad(ggm400{7}/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});

start   = 200;
fin     = 1478;

% Interval & Input function
TimeArray = TimeGPS(start:fin) - TimeGPS(start);
DG1_res        = DG1_res(start:fin); 
DG2_res        = DG2_res(start:fin); 
DG3_res        = DG3_res(start:fin); 

DG1_low        = DG1_low(start:fin); 
DG2_low        = DG2_low(start:fin); 
DG3_low        = DG3_low(start:fin); 

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
    
    DG_res = [DG1_res; DG2_res; DG3_res];
    DG_low = [DG1_low; DG2_low; DG3_low];
else
    DG_res = DG3_res;
    DG_low = DG3_low;
end

% LS estimation
WCoeff_res = Psi \ DG_res;
WCoeff_low = Psi \ DG_low;

disp(['Nr. of wavelet-coef: ',num2str(max(size(WCoeff_res)))])

disp(['Analysing SVD...'])
s=svd(Psi);
if s(1)/s(end)>10^16
disp(['Warning: matrix is ill-conditioned'])    
end
disp(['Condition nr: ',num2str(s(1)/s(end))])


% Wavelet-reconstruction
DG_est_res = Psi * WCoeff_res;

DG1_est_res = DG_est_res(1:end/3, :);
DG2_est_res = DG_est_res(end/3+1:2*end/3, :);
DG3_est_res = DG_est_res(2*end/3+1:end, :);
% DG_est3 = DG_est;

DG_est_low = Psi * WCoeff_low;

DG1_est_low = DG_est_low(1:end/3, :);
DG2_est_low = DG_est_low(end/3+1:2*end/3, :);
DG3_est_low = DG_est_low(2*end/3+1:end, :);
% DG_est3 = DG_est;

figure(1)
subplot(2, 1, 1)
hold on;
plot(TimeArray,DG1_res)
plot(TimeArray, DG1_est_res,'r')
title(['dg1 residium reconstruction. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
legend('true','estimate')
grid on;
xlabel('Time(s)')

subplot(2, 1, 2)
hold on;
plot(TimeArray,DG1_low)
plot(TimeArray, DG1_est_low)
legend('true','estimate')
grid on;
title(['dg1 low reconstruction. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
xlabel('Time(s)')

figure(2)
subplot(2, 1, 1)
hold on;
plot(TimeArray,DG2_res)
plot(TimeArray, DG2_est_res,'r')
title(['dg1 residium reconstruction. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
legend('true','estimate')
grid on;
xlabel('Time(s)')

subplot(2, 1, 2)
hold on;
plot(TimeArray,DG2_low)
plot(TimeArray, DG2_est_low)
title('\dg1 low')
legend('true','estimate')
grid on;
title(['dg1 low reconstruction. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
xlabel('Time(s)')


figure(3)
subplot(2, 1, 1)
hold on;
plot(TimeArray,DG3_res)
plot(TimeArray, DG3_est_res,'r')
title(['dg1 residium reconstruction. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
legend('true','estimate')
grid on;
xlabel('Time(s)')

subplot(2, 1, 2)
hold on;
plot(TimeArray,DG3_low)
plot(TimeArray, DG3_est_low)
title('\dg1 low')
legend('true','estimate')
grid on;
title(['dg1 low reconstruction. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
xlabel('Time(s)')


% delta = j_max - j_min +1;
% WCoeff1 = WCoeff_res(1:length(WCoeff_res)/3);
% figure(4)
% for i = 1:delta
%     subplot(delta,1,i)
%     WCoeffSub = WCoeff1((i-1)/delta*end+1:i/delta*end);
%     plot(WCoeffSub); hold on;
%     plot(WCoeffSub,'.')
%     title(['\Delta g_1 coeffiсients, coeffiсientslevel=',num2str(j_min+i-1)])
%     xlabel('Number of coefficient')
% end
% 
% WCoeff2 = WCoeff_res(length(WCoeff_res)/3+1:2*length(WCoeff_res)/3);
% figure(5)
% for i = 1:delta
%     subplot(delta,1,i)
%     WCoeffSub = WCoeff2((i-1)/delta*end+1:i/delta*end);
%     plot(WCoeffSub); hold on;
%     plot(WCoeffSub,'.')
%     title(['\Delta g_2 coeffiсients, level=',num2str(j_min+i-1)])
%     xlabel('Number of coefficient')
% end
% 
% WCoeff3 = WCoeff_res(2*length(WCoeff_res)/3+1:end);
% figure(6)
% for i = 1:delta
%     subplot(delta,1,i)
%     WCoeffSub = WCoeff3((i-1)/delta*end+1:i/delta*end);
%     plot(WCoeffSub); hold on;
%     plot(WCoeffSub,'.')
%     title(['\Delta g_3 coeffiсients, level=',num2str(j_min+i-1)])
%     xlabel('Number of coefficient')
% end
