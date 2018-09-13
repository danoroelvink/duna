function [u10t,u10t0]=vegetation_critical_velocity(D,Cv,k)
% Compute critical velocity as a function of grain size,
% percentage of plant cover (C) and vegetation shape (k)
% Buckley, R. 1987
rho=1.25;
rhos=2650;
g=9.81;
A=0.08;
kappa=0.4;
z0=D/30;
z=10;

% k=0.018 for erected vegetation shape
k=0.046 %for rounded vegetation shape

ustart=A*sqrt(g*D*((rhos-rho)/rho));

u10t0=ustart/kappa.*log(z./z0);

% simple vegetation effect

u10t=u10t0/(1-k*Cv*100);





