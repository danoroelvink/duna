function [u10t0,ustart]=critical_velocity(p)
% Compute critical wind velocity as a function of D50 grain size,
%% General inputs
%[p] = duna_input();

ustart=sqrt(p.g*p.D50*((p.rhos-p.rho)/p.rho))*p.A;

u10t0=ustart/p.kappa.*log(p.z./p.z0);
