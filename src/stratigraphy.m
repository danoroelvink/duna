function [ xtdata,zstrat,time ] = stratigraphy( x,zbt,t )
%STRATIGRAPHY create and plot stratigraphy from
nt=size(zbt,1);
%nt=length(t);
nx=length(x);%size(x,1);
for it=1:nt
    if it>1
        for ix=1:nx
            for itt=it-1:-1:1
                %% Move eroded layers
                if zbt(itt,ix)>zbt(it,ix)
                    zbt(itt,ix)=zbt(it,ix);
                end
            end
        end
    end
end

%% Create matrices with xt and zbt as 'grid' and time as the colored variable
xtdata=zeros(size(zbt));
time=zeros(size(zbt));

for it=1:nt
    xtdata(it,:)=x;
    time(it,:)=t(it);
end
zstrat=zbt;
end

