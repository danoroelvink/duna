function [p]=duna_input();
% function to input all parameter values for duna
%% duna parameters
p.plotinterval = 1000000;       % Plot interval
p.xplotmin = -450;
p.xplotmax = 0;
p.D50 = 0.00035;                % Grain size
p.rhos = 2650;                  % density of sand
p.por = 0.4;                    % porosity of sand
p.p1 = 0.9 ;%0.95               % fraction of sand susceptible to aeolian transport
p.dzmix = 0.05;                 % thickness of layer for plastering effect
p.nyear=1;                      % number of years to simulate
p.morfac=1;                     % morphological factor
p.dt=7200;                      % default time step
p.zbmin_sup=0;                  % min bed level for supply
p.zbmax_sup=3.5;                % max bed level for supply
p.supply=0.0;                   % supply in m/day
p.M2amp=1.5;                    % M2 tidal amplitude
p.M2freq=2*pi/740/60;           % M2 tidal frequency
p.runup=2;                      % runup 
p.Tm=4*3600;                    % time scale of dropping of moisture content
p.Tg=90*24*3600;                % time scale of vegetation growth
p.T=1;                          % time scale of adaptation of concentration or pick up
p.Cb=2;                         % constant in formulation for equilibrium sediment concentration
p.leeslope=1/10;                % slope of lee zone behind dune
p.duran=16;                     % dimensionless “roughness factor,” 
                                % effectiveness of vegetation in slowing down flow (Duran and Moore, 2013)
%% Kroy model
p.alfa=3;                       % 0.45 times typical length/heigth ratio
p.beta=1;                       % parameter beta in Kroy model
%% Critical velocity
p.rho=1.25;                     % density of air
p.g=9.81;                       % acceleration of gravity
p.A=0.08;                       % calibration factor; A = 1 ???
p.kappa=0.4;                    % von Karman constant
p.z=10;                         % height where wind speed is given
p.ks=0.015;
p.z0=p.ks/30;
p.k=0.9;                        % constant parameter; for small erect plants k=0.018, for small rounded stemless plants k=0.046
p.tanalpha=tan(35*pi/180);      % tangent of angle of repose

%% Vegetation
%p.wind_reduction='duran';
p.wind_reduction='buckley';
p.zbcrit=4.0;
p.veginit = 1
p.zbvegtab=[-1000    4.0   5.2    7  1000];   
p.Hvvegtab=[    0    0.0   0.3    .6    1];
p.cvvegtab=[    0    0.0   0.3    .6   .9]; 

