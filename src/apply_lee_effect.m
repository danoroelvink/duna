function [ Cu ] = apply_lee_effect(x,zb,slope,Cu,phi);
%% APPLY_LEE_EFFECT Set Cu to 0 in lee zones
slope_cor=slope./(max(abs(cos(phi)),1.e-6));
sep=zb;
if cos(phi(1))>=0
    for i=2:length(x)
        sep(i)=sep(i-1)-slope_cor(i)*(x(i)-x(i-1));
        if sep(i)>zb(i);
            Cu(i)=0;
        else
            sep(i)=zb(i);
        end
    end
else
    for i=length(x)-1:-1:1
        sep(i)=sep(i+1)-slope_cor(i)*(x(i+1)-x(i));
        if sep(i)>zb(i);
            Cu(i)=0;
        else
            sep(i)=zb(i);
        end
    end
end
