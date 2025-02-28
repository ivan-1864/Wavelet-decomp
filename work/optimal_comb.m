Filename_ggm = '../data/XGM2019_1000.dat';
fid        =  fopen(Filename_ggm,'r');              
ggm1000        =  textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',44);  
fclose(fid);

Filename_ggm = '../data/XGM2019_400.dat';
fid        =  fopen(Filename_ggm,'r');              
ggm400        =  textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',44);  
fclose(fid);

TimeGPS = ggm400{1};


DG3_high = 10^(-5)*(ggm1000{5}-ggm400{5});
DG1_high = -deg2rad((ggm1000{6}-ggm400{6})/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});
DG2_high = -deg2rad((ggm1000{7}-ggm400{7})/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});

DG3_low = 10^(-5)*(ggm400{5});
DG1_low = -deg2rad(ggm400{6}/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});
DG2_low = -deg2rad(ggm400{7}/ 3600) .* Geodesy_NormalGravity(ggm400{3}, ggm400{4});

% Interval & Input function
TimeArray = TimeGPS(start:fin) - TimeGPS(start);
DG1_high        = DG1_high(start:fin); 
DG2_high        = DG2_high(start:fin); 
DG3_high        = DG3_high(start:fin); 

DG1_low        = DG1_low(start:fin); 
DG2_low        = DG2_low(start:fin); 
DG3_low        = DG3_low(start:fin); 

Time_fin  = TimeArray(end);


df = zeros(0, 5);

for j_min = -20:1:0
    disp(num2str(j_min))
    for delta = 0:20
        for k_max = 20:20:400
            if 3*(delta+1) * k_max <= 600
                j_max = j_min + delta;
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
            
            DG_high = [DG1_high; DG2_high; DG3_high];
            DG_low = [DG1_low; DG2_low; DG3_low];

            % LS estimation
            WCoeff_high = Psi \ DG_high;
            WCoeff_low = Psi \ DG_low;

            n = length(WCoeff_high);
            
            % Wavelet-reconstruction
            DG_est_high = Psi * WCoeff_high;
            DG_est_low = Psi * WCoeff_low;

            rms_high = rmse(DG_high, DG_est_high);
            rms_low = rmse(DG_low, DG_est_low);
            
            df = [df; j_min, j_max, k_max, n, rms_high, rms_low];
            end
        end
    end
end

header = {'j_min', 'j_max', 'k_max', 'n', 'rms_high', 'rms_low'};
ds = dataset({df,header{:}});

% ds_fil = ds(ds.rms < 10^(-6), header);
% ds_fil = sortrows(ds_fil,'n','ascend');

% для 400
% -11	-8	20 -- оптимально по rmse
% -9	-9	100 -- интересно по уровням
% -12	-7	20 -- наибольший размах

% для 1000
% -9	-9  100 -- оптимально по rmse
% -10	-8  40 -- несколько уровней
