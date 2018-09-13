function [xrunup]=getrunup(fname)
%% get the x position of maximum runup, in the original x coordinate
%% for the case of a 1D XBeach run

x=nc_varget(fname,'globalx');
y=nc_varget(fname,'globaly');
alfa=atan2d(y(end)-y(1),x(end)-x(1));
x=x*cosd(alfa)+y*sind(alfa);

zs=nc_varget(fname,'zs');
zb=nc_varget(fname,'zb');
h=zs-zb;
wet=max(h)>0.01;
irunup=find(wet,1,'last');
xrunup=x(irunup);
