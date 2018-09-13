function [xgrid zgrid] = makeprofile(xdata_all,zdata_all,zbmin,zbmax,dxmin,dxmax,phibeach);
hmin=0.1;

%% look up points where zbmin and zbmax are crossed for the first time
for i=1:length(xdata_all)
    if zdata_all(i)>=zbmin
        imin=i-1;
        break
    end
end
imin=max(imin,1);
imax=length(xdata_all);
for i=imin:length(xdata_all)
    if zdata_all(i)>=zbmax
        imax=i;
        break
    end
end

xdata=xdata_all(imin:imax);
zdata=zdata_all(imin:imax);
xgrid(1)=xdata(1);
zgrid(1)=zdata(1);
h0=-zdata(1);
xgrid(2)=xgrid(1)+dxmax;
i=1;
dx=dxmax;
while xgrid(i)+dx<xdata(end)
    i=i+1;
    h=max(-zgrid(i-1),hmin);
    dx=dxmax*sqrt(h/h0);
    dx=max(dx, dxmin);
    xgrid(i)=xgrid(i-1)+dx;
    zgrid(i)=interp1(xdata,zdata,xgrid(i),'linear','extrap');
    
end
%zgrid(1:3)=zgrid(3);
% xgrid=[xdata(1):dx:xdata(end)];
% zgrid=interp1(xdata,zdata,xgrid);
alfa=mod(270-phibeach,360);

fi=fopen('griddata.txt','wt');
fprintf(fi,'nx      = %3i \n',length(xgrid)-1);
fprintf(fi,'ny      = %3i \n',0);
fprintf(fi,'vardx   = 1 \n',1);
fprintf(fi,'dy      = %6.1f \n',100);
fprintf(fi,'xori    = %10.2f \n',0);
fprintf(fi,'yori    = %10.2f \n',0);
fprintf(fi,'alfa    = %10.2f \n',alfa);
fprintf(fi,'xfile   = x.dep \n');
fprintf(fi,'depfile = profile.dep \n');
fprintf(fi,'posdwn  = -1 \n');
fclose(fi);

fi=fopen('x.dep','wt');
for i=1:1
    fprintf(fi,'%7.3f ',xgrid);
    fprintf(fi,'\n');
end
fclose(fi);
fi=fopen('profile.dep','wt');
for i=1:1
    fprintf(fi,'%7.3f ',zgrid);
    fprintf(fi,'\n');
end
fclose(fi);

if 0
    figure(3);
    plot(xdata,zdata,'-',xgrid,zgrid,'o','linewidth',2);
    legend('data','grid');
end
