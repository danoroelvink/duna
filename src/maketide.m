function [zs0] = maketide(MLWS,MHWS,tmax)
%MAKETIDE genrate a sinusoidal tide time series with a random phase

t=[0:600:tmax+1200];
omega=2*pi/(12.5*3600);
phase=rand*2*pi;
ztide=.5*(MHWS-MLWS)*sin(omega*t-phase);
clear out;
out(:,1)=t;
out(:,2)=ztide;
fi=fopen('tide.txt','w');
fprintf(fi,'%8.0f %8.2f\n',out');
fclose(fi);
zs0=ztide(1);
end

