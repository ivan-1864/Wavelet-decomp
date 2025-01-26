DataFile = '../Data/Phase_L1.VEL';
TimeCol = 1;
DataCol = 4;

% number of wavelets on freq-axis     
j_max = 4;

% number of wavelets on t-axis
k_max = 25;

% -----------end cfg-----------------------------

D = importdata(DataFile);

Time = D.data(:, TimeCol);
Hei = D.data(:, DataCol);

clear D

% mask = find(data(:, 1) >= 16130 & data(:, 1) <= 18549);
% data = data(mask, :);


j_k_max = k_max * 2^j_max;

Time_st = min(Time);
Time_fin = max(Time);




% normalize time
time_wev = (Time - Time_st) / (Time_fin - Time_st) * j_k_max;

Psi = zeros(length(Time), k_max*(2^(j_max+1) - 1));
for t = 1 : length(Time)
    shift = 1;
    time = time_wev(t);
    for j = 0 : j_max
        for k = 0 : k_max * 2^(j_max-j)
          Psi(t, k+shift) = wavel_trf(j, k, time);
        end
        shift = shift + k_max * 2^(j_max-j);
    end
end

d = Psi \ Delta_g  ;

Delta_g_app = Psi * d;

% shift = 1;
% for j = 0 : j_max
%     k = 0 + shift : k_max * 2^(j_max-j) + shift;
%     Delta_g_app_lev() = Psi(:, k) * d(k);
%     figure(j+1)
%     hold
%     plot(Time, Delta_g_app_lev)
%     shift = shift + k_max * 2^(j_max-j);
% end


% for j = 0 : j_max
%     ind = j * (k_max + 1) +1 : (j + 1)  * (k_max + 1); 
%     Delta_g_part = Psi(:, ind)*d(ind);
%     figure(1)
%     hold
%     plot(Time, Delta_g_part)
% end
ind = 401:776;
Delta_g_app_1 = Psi(:, ind)*d(ind);

figure(10)
hold
plot(Time, Delta_g_app)
plot(Time, Delta_g_app_1)
plot(Time, Delta_g)