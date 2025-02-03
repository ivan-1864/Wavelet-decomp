% GNSS solutions 
%
% Vadim Vyazmin (NavLab, MSU)
% 09/11/2020
%----------------------------------------------------------------------------
clc; close all; 
M=12;

FileName_GPS = '../Data/Phase_L1.VEL';


%--------------- SR2NAV solution: ------------------
D            = importdata(FileName_GPS); % Time[s] Lat[rad] Lon[rad] Hei[m] RmsPos V_E V_N V_UP RmsVel SVs Type
disp(['Phase-vel file:  ',FileName_GPS])
gnss         = D.data;
Dtxt = D.textdata;
Date = '2024:01:30';%Dtxt{1}(26:35);


%------------------------------------------------------------------------------------
i0_gps=1;
in_gps=max(size(gnss));

dt_gps    = mean(diff(gnss(1:100,1)));
TimeGPS   = gnss(:,1);
Lat_gps     = rad2deg(gnss(i0_gps:in_gps,2));  % deg
Lon_gps     = rad2deg(gnss(i0_gps:in_gps,3));  % deg
Hei_gps     = gnss(i0_gps:in_gps,4);           % m
Ve_gps      = gnss(i0_gps:in_gps,6);           % m/s
Vn_gps      = gnss(i0_gps:in_gps,7);           % m/s
Vup_gps     = gnss(i0_gps:in_gps,8);           % m/s
RmsPos      = gnss(i0_gps:in_gps,5);           % m
RmsVel      = gnss(i0_gps:in_gps,9);           % m/s
SVs         = gnss(i0_gps:in_gps,10); 
Type_sln    = gnss(i0_gps:in_gps,11);          %
clear gnss

figure('Name','Sr2Nav - Hei'); clf;
plot(TimeGPS,Hei_gps)
h = title(['Height from GPS. Flight ', num2str(Date)]);
h1 = xlabel('Time (s)');
h2 = ylabel('Height above ell. (m)');
grid on


figure('Name','Sr2Nav - Traj'); clf;
plot(Lon_gps,Lat_gps)
h = title(['Trajectory from GPS. Flight ', num2str(Date)]);
h1 = xlabel('Lon (deg)');
h2 = ylabel('Lat (deg)');
grid on


figure('Name','Sr2Nav - QC'); clf;
plot(TimeGPS,RmsPos)
hold on
plot(TimeGPS,RmsVel,'r')
hold on
plot(TimeGPS,Type_sln,'g')
hold on
plot(TimeGPS,Type_sln,'g.')
hold on
plot(TimeGPS,RmsVel,'r.')
h = title(['GPS solution QC (Diff mode). Flight ', num2str(Date)]);
h1 = xlabel('Time (s)');
legend('RmsPos','RmsVel','Type')
grid on

figure('Name','Sr2Nav - Vel'); clf;
T_const = 0; %s  Days of week
plot(TimeGPS+T_const,[Ve_gps,Vn_gps,Vup_gps])
title(['Velocity from GPS. Flight ', num2str(Date)])
ylabel('(m/s)')
xlabel('Time (s)');
legend('V_E','V_N','V_{Up}')
grid on


figure('Name','Sr2Nav - Horiz_vel'); clf;
plot(TimeGPS,sqrt(Ve_gps.^2+Vn_gps.^2))
title(['Horizontal velocity from GPS. Flight ', num2str(Date)])
ylabel('(m/s)')
xlabel('Time (s)');
grid on


figure('Name','Sr2Nav - SVs'); clf;
plot(TimeGPS, SVs,'k')
title(['SVs (Diff mode). Flight ', num2str(Date)]);
xlabel('Time (s)');
grid on


% % Trajectory in Gauss-Krueger
% Rgk = Geodesy_GaussKruger( Lat_gps, Lon_gps, min(Lat_gps), min(Lon_gps) ); 
% 
% figure('Name','Sr2Nav - Traj_GK'); clf;
% plot(Rgk(:,1)/1000,Rgk(:,2)/1000,'k','Linewidth',2.0);
% hold on
% plot(Rgk(1:1,1)/1000,Rgk(1:1,2)/1000,'r*','Linewidth',5.0);
% h = title(['Trajectory from GPS (Gauss-Krueger). Flight ', num2str(Date)]);
% h1 = xlabel('Longitude (km)');
% h2 = ylabel('Latitude (km)');
% set(h ,'fontsize',M,'fontname','Arial Narrow');
% set(h1,'fontsize',M,'fontname','Arial Narrow');
% set(h2,'fontsize',M,'fontname','Arial Narrow');
% set(gca,'fontsize',M,'fontname','Arial Narrow'); 
% box on
% grid on
% %axis equal