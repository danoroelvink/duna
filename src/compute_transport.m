function [ Cc,Sx,Cu ] = compute_transport( Cc,x,dx,dt,u10,u10x,u10t0,dir,Cu,Cb,T,p1,wet,OffshoreBound,projectpr)
%% Compute sand transport by solving 1D advection equation

switch projectpr
    case 'yes'
% % SOS projecting x to the wind direction
factor=1/cos(dir);
    otherwise
factor=1.0;
end
Sx=zeros(size(Cc));

% Critical fetch length from Delgado-Fernandez(2010)
Fc=max(0,4.38*u10-8.23); 
if u10x(1)>0  % onshore wind
    switch OffshoreBound
        case 'closed'
    Cc(1)=0.0;%changing boundary concentration to zero
        case 'open'    
    Cc(1)=Cu(1); 
    end
    Sx(1)=Cc(1)*u10x(1);
    xx=0;
    io=find(wet==1,1,'last');
    if io>2
        Sx(2:io)=0.0;
        Cc(2:io)=0.0;  
    else
        io=2;
    end
    
%    sin(pi/2*F/Fc)
    for i=io:length(Cc) 
        xx=xx+factor*abs(x(i-1)-x(i));
        if xx<Fc
        Cfetch=sin(pi/2*xx/Fc); % from Delgado-Fernandez(2010)
        else
        Cfetch=1.0;    
        end
        Cu(i)=Cu(i)*Cfetch;

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