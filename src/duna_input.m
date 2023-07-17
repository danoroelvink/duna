function [p]=duna_input(storm,OffshoreBound)
% function to input all parameter values for duna
%% duna parameters
p.projectprofile='yes';
p.plotinterval = 1000000;              % Plot interval
p.xplotmin = -450;
p.xplotmax = 0;
p.D50 = 0.00035;                 % Grain size
p.rhos = 2650;                  % density of sand
p.por = 0.4;                    % porosity of sand

switch num2str(storm)
    case '1'
p.runup=4;                    % runup during Emma storm
    otherwise
p.runup=0.75;                   % typical runup   
end
p.p1 = 0.9;             % fraction of sand susceptible to aeolian transport 
p.dzmix = 0.05;                 % thickness of layer for plastering effect
p.nyear=1;                      % number of years to simulate
p.morfac=1;                     % morphological factor
p.dt=300;                      % default time step
switch OffshoreBound
case 'open'
p.zbmin_sup=-1;                  % min bed level for supply
p.zbmax_sup=0.2;                % max bed level for supply    
p.supply=0.0;                   % supply in m/day
case 'closed'
p.zbmin_sup=100.0;             % min bed level for supply
p.zbmax_sup=99.0;              % max bed level for supply
p.supply=0.1;                   % supply in m/day
end

p.M2amp=1.5;                    % M2 tidal amplitude
p.M2freq=2*pi/740/60;           % M2 tidal frequency
% p.vegtoler=0.2;                 % vegetation tolerance to burial (m)
p.Tm=4*3600;                    % time scale of dropping of moisture content
p.Tg=360*24*3600/2.;                % time scale of vegetation growth
p.T=0.5;                         % time scale of adaptation of concentration or pick up
p.Cb=2;                         % constant in formulation for equilibrium sediment concentration
p.leeslope=4/10.;               % slope of lee zone behind dune
p.duran=16;                     % dimensionless “roughness factor,” (was 16)
                                % effectiveness of vegetation in slowing down flow (Duran and Moore, 2013)
%% Kroy model
% p.alfa=3;                     % 0.45 times typical length/heigth ratio
% p.beta=1;                     % parameter beta in Kroy model
p.alfa=3;                       %   --> sensitivity to slope changes
p.beta=0.2;                     %   --> balance between stoss/lee : high/low
p.alfa2=2;                      %   --> sensitivity to slope changes
p.beta2=1;                      %   --> balance between stoss/lee : high/low
%% Critical velocity
p.rho=1.25;                     % density of air
p.g=9.81;                       % acceleration of gravity
p.A=0.08;                       % calibration factor; A = 1 ???
p.kappa=0.4;                    % von Karman constant
p.z=10;                         % height where wind speed is given
p.ks=0.015;
p.z0=p.ks/30;
p.k=0.2;                        % constant parameter; for small erect plants k=0.018, for small rounded stemless plants k=0.046 %we were using 0.9
p.k2=0.45;                      % buckley constant for crest-lee
p.k3=0.8;                       % buckley constant for lee
p.tanalpha=tan(35*pi/180);      % tangent of angle of repose

%% Vegetation
% p.wind_reduction='duran';
p.wind_reduction='buckley';
p.veginit = 1;
p.zbcrit=5.5; 
p.zbvegtab=[-1000    5.4   5.5      6.       6.5     7.      1000];   
p.Hvvegtab=[    0    0.0   0.15     0.15     0.2     0.3     1];
p.cvvegtab=[    0    0.0   0.20     0.40	 0.7     0.7     1]; 
