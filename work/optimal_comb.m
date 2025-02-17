read_ggm;

start = 200;
fin = 1478;
TimeArray = TimeGPS(start:fin) - TimeGPS(start);
DG        = DG(start:fin); 
Time_fin  = TimeArray(end);

df = zeros(0, 5);

for j_min = -20:1:0
    disp(num2str(j_min))
    for delta = 0:20
        for k_max = 20:20:400
            if delta * k_max <= 500
                j_max = j_min + delta;
                main;
            end
        end
    end
end