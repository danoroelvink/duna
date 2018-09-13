function [zs0] = maketide2(MLWS,MHWS,tmax)
%MAKETIDE genrate a sinusoidal tide time series 
T=12.5*3600;
ntides=round(5*tmax/T)
omega=2*pi/tmax*ntides;
t=[0:60:tmax+660];
phase=600*omega;
ztide=-.5*(MHWS-MLWS)*cos(omega*t-phase)+.5*(MLWS+MHWS);
clear out;
out(:,1)=t;
out(:,2)=ztide;
fi=fopen('tide.txt','w');
fprintf(fi,'%8.0f %8.2f\n',out');
fclose(fi);
zs0=ztide(1);
% figure;
% plot(t/3600,ztide);
end

