function [ tau_over_tau0 ] = kroy( x,zb,alfa,beta )
% Kroy computes shear stress/undisturbed shear stress ratio over arbitrary 
% dune profile zb given on equidistant grid x
% ref. Kroy,Sauermann,Hermann: Minimal Model for Sand Dunes, eq. 4

%% Initialize
dzbdx=zeros(size(x));
tau_over_tau0=zeros(size(x));

%% Compute bed gradients
dx=x(2)-x(1);
dzbdx(2:end-1)=(zb(3:end)-zb(1:end-2))/2/dx;

%% Compute convolution integral
nx=length(x);
for i=1:nx
    integ=0;
    for j=i-nx:i-1
        if j~=0
            integ=integ+dzbdx(i-j)/(j*pi);
        end
    end
    tau_over_tau0(i)=alfa*(integ+beta*dzbdx(i))+1;
end
tau_over_tau0=max(tau_over_tau0,0.1);
end

