% Wavelet-approximation by LSE
%
% Version: one-level. Date: 31.01.2025
% close all; clc;

% % Wavelet settings
% j_min = -8;
% j_max = -7;       % resolution level 
% k_max = 200;      % nr. of wavelets on t-axis

% GNSS data
% read_gnss_sr2nav; 
% read_ggm;

% close all

% Interval & Input function
% start = 200;
% fin = 1478;
% TimeArray = TimeGPS(start:end) - TimeGPS(start);
% DG        = DG(start:end); 
% Time_fin  = TimeArray(end);


% Matrix of wavelet values

Psi = zeros( length(TimeArray), (j_max - j_min + 1) * k_max );
for j = j_min : j_max
    for k = 1 : k_max
       dT  = (2^j)*Time_fin/k_max;  % step of wavelet grid
       t_k = k * dT;                % knot in grid
       Psi(:,(j - j_min) * k_max + k) = wavel_trf(j,t_k,TimeArray);   
    end
end

% LS estimation
WCoeff = Psi \ DG;
% disp(['Nr. of wavelet-coef: ',num2str(max(size(WCoeff)))])

% disp(['Analysing SVD...'])
% s=svd(Psi);
% if s(1)/s(end)>10^16
% disp(['Warning: matrix is ill-conditioned'])    
% end
% disp(['Condition nr: ',num2str(s(1)/s(end))])
n = length(WCoeff);

% Wavelet-reconstruction
DG_est = Psi * WCoeff;
rms = rmse(DG, DG_est);

df = [df; j_min, j_max, k_max, n, rms];


% figure(1)
% plot(WCoeff); hold on;
% plot(WCoeff,'.')
% grid on;
% title(['Estimated wavelet-coefficients. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
% xlabel('Number of coefficient')
% 
% figure(2)
% plot(TimeArray,DG)
% hold on;
% plot(TimeArray, DG_est,'r')
% legend('true','estimate')
% grid on;
% title(['Reconstructed input function. MHat wavelet, level=',num2str(j_min),'-',num2str(j_max)])
% xlabel('Time(s)')
% 
% disp(['RMSE = ', num2str(round(rmse(DG, DG_est), 3))])

