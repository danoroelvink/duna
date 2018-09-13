function [ Cu ] = compute_Cu( umagx,u10t0 )
% Function to compute equilibriunm concentration
Cu=1.5e-4*(max(umagx-u10t0,0)).^3./umagx;

end

