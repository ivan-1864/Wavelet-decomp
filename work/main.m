% Wavelet-approximation by LSE
%
% Version: one-level. Date: 31.01.2025
close all; clc;

% Wavelet settings
j_min = -12;
j_max = -7;       % resolution level 
k_max = 20;      % nr. of wavelets on t-axis

% GNSS data
% read_gnss_sr2nav; 
% close all

Filename_ggm = '../data/XGM2019_400.dat';
fid        =  fopen(Filename_ggm,'r');              
ggm        =  textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',44);  
fclose(fid);

TimeGPS = ggm{1};

Lon = ggm{2};
Lat = ggm{3};

DG3 = 10^(-5)*ggm{5};
DG1 = -deg2rad(ggm{6}/ 3600) .* Geodesy_NormalGravity(ggm{3}, ggm{4});
DG2 = -deg2rad(ggm{7}/ 3600) .* Geodesy_NormalGravity(ggm{3}, ggm{4});

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


figure(1)
plot(TimeArray,(DG1-DG_est1)*10^5)
hold on;
% plot(TimeArray, DG_est1,'r')
% legend('true','estimate')
grid on;
title(['Reconstructed input function. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
xlabel('Time(s)')

figure(2)
plot(TimeArray,(DG2-DG_est2)*10^5)
hold on;
% plot(TimeArray, DG_est2,'r')
% legend('true','estimate')
grid on;
title(['Reconstructed input function. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
xlabel('Time(s)')

figure(3)
plot(TimeArray,(DG3-DG_est3)*10^5)
hold on;
% plot(TimeArray, DG_est3,'r')
% legend('true','estimate')
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
    title(['Estimated wavelet-coefficients. MHat wavelet, level=',num2str(j_min+i-1)])
    xlabel('Number of coefficient')
end

WCoeff2 = WCoeff(length(WCoeff)/3+1:2*length(WCoeff)/3);
figure(5)
for i = 1:delta
    subplot(delta,1,i)
    WCoeffSub = WCoeff2((i-1)/delta*end+1:i/delta*end);
    plot(WCoeffSub); hold on;
    plot(WCoeffSub,'.')
    title(['Estimated wavelet-coefficients. MHat wavelet, level=',num2str(j_min+i-1)])
    xlabel('Number of coefficient')
end

WCoeff3 = WCoeff(2*length(WCoeff)/3+1:end);
figure(6)
for i = 1:delta
    subplot(delta,1,i)
    WCoeffSub = WCoeff3((i-1)/delta*end+1:i/delta*end);
    plot(WCoeffSub); hold on;
    plot(WCoeffSub,'.')
    title(['Estimated wavelet-coefficients. MHat wavelet, level=',num2str(j_min+i-1)])
    xlabel('Number of coefficient')
end
