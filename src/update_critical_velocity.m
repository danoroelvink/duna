function [u10t0]=update_critical_velocity(m,x,zb,ustart,phi,p)
% Updates critical wind velocity as a function of slope angle, moisture
% content and vegetation

%% slope effect (Hardisty and Whitehouse, 1988; de Vries et al., 2012); 
%% slope has to be determined between i and its downwind point, so
%% depending on wind direction
nx=length(x);
if cos(phi(1))>0
    theta(nx)=atan((zb(nx)-zb(nx-1))/(x(nx)-x(nx-1))); % beach slope
    for i=1:nx-1
        theta(i)=atan((zb(i+1)-zb(i))/(x(i+1)-x(i)));
    end
else
    theta(1)=atan((zb(2)-zb(1))/(x(2)-x(1))); % beach slope
    for i=2:nx
        theta(i)=atan((zb(i)-zb(i-1))/(x(i)-x(i-1)));
    end
end

ustart=ustart.*sqrt(max((cos(theta)+cos(phi).*sin(theta)./p.tanalpha),0));
save('critvel0.mat','theta','phi','ustart')
%% moisture content 
% uses the geotechnical mass content (ratio of dry mass)
% m=(moisture*rhow)/(rhos(1 - porosity))

moist = m > 0.0032;      % minimum water content in nature
    ustart(moist)=ustart(moist)+7.5*m(moist); % Hotta et al. 1984
moister = m > 0.064;    % corresponds to the 10% of volume water content (Delgado-Fernandez, 2010)
    ustart(moister)=inf;% set to infinity if the content is above 10%
    
u10t0=ustart/p.kappa.*log(p.z./p.z0);
save('critvel.mat','ustart','u10t0','moist','moister')
