function [zb] = avalan(x,zb,p)
dza=zeros(size(x));
dx=x(2)-x(1);
for i=1:length(x)-1
    dz=zb(i+1)-zb(i);
    if abs(dz)>p.tanalpha*dx
        ddz=.5*(abs(dz)-p.tanalpha*dx)*sign(dz);
        dza(i)=dza(i)+ddz;
        dza(i+1)=dza(i+1)-ddz;
    end
end
zb=zb+dza;

end

