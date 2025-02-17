Filename_ggm = '../data/XGM2019_400.dat';
fid        =  fopen(Filename_ggm,'r');              
ggm        =  textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',44);  
fclose(fid);

start = 200;
fin = 1478;

TimeGPS = ggm{1};
DG = ggm{5};