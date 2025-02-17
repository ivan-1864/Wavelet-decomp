Filename_ggm = '../data/XGM2019_1000.dat';
fid        =  fopen(Filename_ggm,'r');              
ggm        =  textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',44);  
fclose(fid);

TimeGPS = ggm{1};


DG3 = 10^(-5)*ggm{5};
DG1 = -deg2rad(ggm{6}/ 3600) .* Geodesy_NormalGravity(ggm{3}, ggm{4});
DG2 = -deg2rad(ggm{7}/ 3600) .* Geodesy_NormalGravity(ggm{3}, ggm{4});

% Interval & Input function
TimeArray = TimeGPS(start:fin) - TimeGPS(start);
DG1        = DG1(start:fin); 
DG2        = DG2(start:fin); 
DG3        = DG3(start:fin); 
Time_fin  = TimeArray(end);


df = zeros(0, 5);

for j_min = -20:1:0
    disp(num2str(j_min))
    for delta = 0:20
        for k_max = 20:20:400
            if 3*(delta+1) * k_max <= 600
                j_max = j_min + delta;
                main;
            end
        end
    end
end

header = {'j_min', 'j_max', 'k_max', 'n', 'rms1', 'rms2', 'rms3', 'rms'};
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
