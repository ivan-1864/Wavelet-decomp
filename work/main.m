DataFile = '../Data/Phase_L1.VEL';
TimeCol = 1;
DataCol = 4;

% number of wavelets on freq-axis     
j_max = 2;

% number of wavelets on t-axis
k_max = 20;

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

Delta_g  = Hei;


% normalize time
time_wev = (Time - Time_st) / (Time_fin - Time_st) * j_k_max;

Psi = zeros(length(Time), k_max*(2^(j_max+1) - 1));
for t = 1 : length(Time)
    shift = 1;
    time = time_wev(t);
    for j = 0 : j_max
        for k = 0 : k_max * 2^j
          Psi(t, k+shift) = wavel_trf(j, k, time);
        end
        shift = shift + k_max * 2^j;
    end
end

d = Psi \ Delta_g  ;

Delta_g_app = Psi * d;
Delta_g_lev = zeros(length(Time), j_max+1);
shift = 1;
for j = 0 : j_max
    index = shift:k_max * 2^j+shift;
    Delta_g_lev(:, j+1) = Psi(:, index) * d(index);
    shift = shift + k_max * 2^j;
end


figure(1)
hold
plot(Time, Delta_g)
plot(Time, Delta_g_app)
plot(Time, Delta_g_lev)
legend('true', 'comp', '1', '2', '3')
