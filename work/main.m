
% Matrix of wavelet values

Psi = zeros( length(TimeArray), (j_max - j_min + 1) * k_max );
for j = j_min : j_max
    for k = 1 : k_max
       dT  = (2^j)*Time_fin/k_max;  % step of wavelet grid
       t_k = k * dT;                % knot in grid
       Psi(:,(j - j_min) * k_max + k) = wavel_trf(j,t_k,TimeArray);   
    end
end
Psi = [Psi, zeros(size(Psi)), zeros(size(Psi));
       zeros(size(Psi)),  Psi, zeros(size(Psi));
       zeros(size(Psi)), zeros(size(Psi)), Psi;];

DG = [DG1; DG2; DG3];

% LS estimation
WCoeff = Psi \ DG;

n = length(WCoeff);

% Wavelet-reconstruction
DG_est = Psi * WCoeff;

DG_est1 = DG_est(1:end/3, :);
DG_est2 = DG_est(end/3+1:2*end/3, :);
DG_est3 = DG_est(2*end/3+1:end, :);

rms1 = rmse(DG1, DG_est1);
rms2 = rmse(DG2, DG_est2);
rms3 = rmse(DG3, DG_est3);
rms = rmse(DG, DG_est);

df = [df; j_min, j_max, k_max, n, rms1, rms2, rms3, rms];



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

