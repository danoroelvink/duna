function [ Cc,Sx ] = compute_transport( Cc,dx,dt,u10x,u10t0,Cu,Cb,T,p1 )
%% Compute sand transport by solving 1D advection equation

Sx=zeros(size(Cc));

if u10x(1)>0  % onshore wind
    Cc(1)=Cu(1);
    Sx(1)=Cc(1)*u10x(1);
    for i=2:length(Cc)
        if min(Cu(i),Cb)-Cc(i)>0    % plastering effect if eroding
            Cc(i)=(min(Cu(i),Cb)/T*p1(i)+Cc(i)/dt+u10x(i-1)/dx*Cc(i-1))...
                /(1/dt+u10x(i)/dx+1/T*p1(i));
        else
            Cc(i)=(min(Cu(i),Cb)/T+Cc(i)/dt+u10x(i-1)/dx*Cc(i-1))/(1/dt+u10x(i)/dx+1/T);
        end
        Sx(i)=u10x(i)*Cc(i);
    end
else
    Cc(length(Cc))=Cu(length(Cc));
    Sx(length(Cc))=Cc(length(Cc))*u10x(length(Cc));
    for i=length(Cc)-1:-1:1
        if min(Cu(i),Cb)-Cc(i)>0    % plastering effect if eroding
            Cc(i)=(min(Cu(i),Cb)/T*p1(i)+Cc(i)/dt-u10x(i+1)/dx*Cc(i+1))...
                /(1/dt-u10x(i)/dx+1/T*p1(i));
        else
            Cc(i)=(min(Cu(i),Cb)/T+Cc(i)/dt-u10x(i+1)/dx*Cc(i+1))/(1/dt-u10x(i)/dx+1/T);
        end
        Sx(i)=u10x(i)*Cc(i);
    end
end

