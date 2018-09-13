function [zb,hv,cv]= duna (x,zb0,durat,phi0,u10,zbprev,hv,cv,xrunup,durat_nowind);
% Main routine Duna; computes profile change due to wind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%
% x		uniformly spaced cross-shore distance, increasing shoreward (m)
% zb0           initial profile elevation (m)
% durat         duration of wind event (s)
% phi0          wind direction (rad relative to shore normal)
% u10           undisturbed wind speed at 10m (m/s)
% zbprev        profile when last leaving duna
% hv            actual vegetation height
% cv            actual vegetation cover
% xrunup        x location of maximum runup since last duna call
% durat_nowind  duration of no-wind events since last call (to allow growth only)

%% General inputs
[p] = duna_input();
zb=zb0;
dx=x(2)-x(1);

%% Initialize
Cc=zeros(size(x));
m=zeros(size(x));
theta=zeros(size(x));
dzbdt=zeros(size(x));
dz=zeros(size(x));
p1=zeros(size(x))+p.p1;

%% Vegetation maximum as a function of height
Hvmax=interp1(p.zbvegtab,p.Hvvegtab,zb);       % Maximum plant height =f(zb)
Cvmax=interp1(p.zbvegtab,p.cvvegtab,zb);       % Maximum plant height =f(zb)
%% update Cvmax and Hvmax so they are constant or increasing landward
icrit=find(zb>p.zbcrit,1,'first');
for i=icrit:length(x)
    Hvmax(i)=max(Hvmax(i),Hvmax(i-1));
    Cvmax(i)=max(Cvmax(i),Cvmax(i-1));
end
%% Initial values for cv if first time and p.veginit==1
if isfield(p,'veginit')
    if p.veginit==1 & sum(cv)==0;
        hv=Hvmax;
        cv=Cvmax;
    end
end
%% Vegetation growth during past non-Duna events, with total duration durat_nowind
hv=Hvmax+(hv-Hvmax).*exp(-durat_nowind/p.Tg);
cv=Cvmax+(cv-Cvmax).*exp(-durat_nowind/p.Tg);

%% Adjust vegetation to account for morphology changes since last call
hv=max(hv-(zb-zbprev),0);     % Vegetation height reduced by sanding up
hv(zb-zbprev<-0.2)=0;         % Vegetation killed if erosion >0.2 m
hv(x<=xrunup)=0;              % Vegetation killed by runup

%% Assume plant height uniformly distributed between 0 and Hvmax, variable cv. Then:
cv=max((Hvmax-(zb-zbprev))./Hvmax,0).*cv;     % Vegetation density                     reduced by sanding up
cv=min(cv,Cvmax);
cv(zb-zbprev<-0.2)=0;         % Vegetation killed if erosion >0.2 m
cv(x<=xrunup)=0;              % Vegetation killed by runup

%% Adjust age effect p to account for morphology change (reset events)
%p=max((abs(zb-zbprev)>.2),p);

%% Start time loop
dt0=p.dt;
p.nt=max(round(durat/p.dt/p.morfac),1);
p.dt=durat/p.nt/p.morfac;
for it=1:p.nt
    t=(it-1)*p.dt;
    
    %% Total water level
    zs=p.M2amp*cos(p.M2freq*t)*ones(size(x))+p.runup;
    
    %% Wind speed distribution
    if cos(phi0)>0
        % positive wind direction
        u10x=u10*cos(phi0) * sqrt( kroy(x,max(zb,zs),p.alfa,p.beta));
    else
        % negative wind direction: reverse computation direction
        u10x=u10*cos(phi0) *              ...
            sqrt( kroy(x(end)-x(end:-1:1), ...
            max(zb(end:-1:1),zs(end:-1:1)),p.alfa,p.beta));
        u10x=u10x(end:-1:1);
    end
    v10x=u10*sin(phi0);
    umagx=sqrt(u10x.^2+v10x.^2);
    phi=atan2(v10x,u10x);
    
    %% Threshold velocity
    [u10t0,ustart]=critical_velocity(p);
    
    %% Effect of wind erosion/sedimentation on vegetation height and density
    hv=max(hv-dz,0);                   % Vegetation height reduced by sanding up
    hv=min(hv,Hvmax);                  % For case of erosion limit hv to Hvmax
    cv=max((Hvmax-dz)./Hvmax,0).*cv;   % Vegetation density reduced by sanding up;
    cv=min(cv,Cvmax);                  % For case of erosion limit cv to Cvmax
    %% Vegetation maximum as a function of height
    Hvmax=interp1(p.zbvegtab,p.Hvvegtab,zb);       % Maximum plant height =f(zb)
    Cvmax=interp1(p.zbvegtab,p.cvvegtab,zb);       % Maximum plant height =f(zb)
    %% update Cvmax and Hvmax so they are constant or increasing landward
    icrit=find(zb>p.zbcrit,1,'first');
    for i=icrit:length(x)
        Hvmax(i)=max(Hvmax(i),Hvmax(i-1));
        Cvmax(i)=max(Cvmax(i),Cvmax(i-1));
    end
    %% Vegetation growth
    hv= hv+p.dt*(Hvmax-hv)/p.Tg;% Actual vegetation height growth
    cv= cv+p.dt*(Cvmax-cv)/p.Tg;% Actual vegetation density growth
    
    %rhov=min((hv./max(Hvmax,1e-6)).^2,1);          % Vegetation density
    rhov=(hv./Hvmax.^2);
    
    %% Reduction of velocity due to vegetation
    switch p.wind_reduction
        case 'duran'
            u10x=u10x./((1+p.duran.*rhov).^2);
            v10x=v10x./((1+p.duran.*rhov).^2);
            umagx=sqrt(u10x.^2+v10x.^2);
        case 'buckley'
            u10x=u10x.*(1-min(p.k*cv,0.7));
            v10x=v10x.*(1-min(p.k*cv,0.7));
            umagx=sqrt(u10x.^2+v10x.^2);
    end
    if max(abs(umagx))>u10t0
        
        %% Compute water level and wet cells
        wet=zeros(size(x));
        i=1;
        while zs(i)>zb(i)&i<length(zb)
            wet(i)=1;
            i=i+1;
        end
        
        %% Update moisture content
        %m=m-m*dt/Tm;
        m(wet==1)=0.25; % maximum water mass content (porosity*rhow/rhos(1-porosity))
        p1(wet==1)=p.p1; % when wet reset to original fraction
        m=m-m/p.Tm*p.dt;
        
        %% Compute critical velocity
        [u10t0]=update_critical_velocity(m,x,zb,ustart,phi,p);
        
        %% Compute equilibrium concentrations
        [ Cu ] = compute_Cu( umagx,u10t0 );
        
        %% Set Cu to 0 in lee zones
        [ Cu ] = apply_lee_effect(x,zb,p.leeslope,Cu,phi);
        
        if isfield(p,'z_nonero')
            p1(zb<p.z_nonero)=0;
            p1(zb>=p.z_nonero)=1;
        end
        %% Update concentrations
        [ Cc,Sx ] = compute_transport( Cc,dx,p.dt,u10x,u10t0,Cu,p.Cb,p.T,p1 );
        
        %% Update bed level
        if cos(phi0)>0
            for i=2:length(x)
                dzbdt(i)=-(Sx(i)-Sx(i-1))/dx/(1-p.por)/p.rhos;
            end
        else
            for i=1:length(x)-1
                dzbdt(i)=-(Sx(i+1)-Sx(i))/dx/(1-p.por)/p.rhos;
            end
        end
        % add a supply
        for i=1:length(x)
            if zb(i)>p.zbmin_sup&zb(i)<p.zbmax_sup
                im1=max(i-1,1);
                ip1=min(i+1,length(x));
                if p.supply>0
                    slope=(zb(ip1)-zb(i))/(x(ip1)-x(i));
                else
                    slope=(zb(i)-zb(im1))/(x(i)-x(im1));
                end
                %dzbdt(i)=dzbdt(i)+p.supply/24/3600*slope;  % supply in m/day
                dzbdt(i)=p.supply/24/3600*slope;  % supply in m/day
            end
            if zb(i)>p.zbmax_sup
                break
            end
        end
        dzbdt(1)=0;
        
        dz=dzbdt*p.morfac*p.dt;
        zb=zb+dz;
        
        %% Avalanching
        zbefore=zb;
        zb=avalan(x,zb,p);
        dz=dz+zb-zbefore;
        
        %% Plastering effect: p1 is fraction of sand susceptible to aeolian transport
        for i=1:length(p1)
            if dz(i)<0
                p1(i)=max(p1(i)+dz(i)/p.dzmix.*(1-p.p1),0);
            else
                p1(i)=max(p1(i)+dz(i)/p.dzmix.*(1-p1(i)),0);
            end
        end
        
        %% Optional plotting
        if mod(it,p.plotinterval)==0
            p.xplotmin=x(1);
            p.xplotmax=x(end);
            subplot(611)
            xgrass=[x;x];
            grass=[zb;zb+hv];
            plot(xgrass,grass,xgrass',grass','color',[0 .5 0],'linewidth',2);
            hold on
            zs(~wet)=nan;
            plot(x,zs,'b',x,zb0,'k',x,zb,'r','linewidth',2);title('zs,zb0,zb,veg')
            hold off
            set(gca,'xlim',[p.xplotmin p.xplotmax],'ylim',[-1 10])
            
            subplot(612)
            plot(x,m,'b',x,p1,'r','linewidth',2);title('m, p1')
            set(gca,'xlim',[p.xplotmin p.xplotmax])
            
            subplot(613);
            plot(x,umagx,'k',x,abs(u10x),'b',x,u10t0,'r','linewidth',2);title('umag,u10,u10t0')
            set(gca,'xlim',[p.xplotmin p.xplotmax])
            set(gca,'ylim',[-max(umagx) max(umagx)])
            
            subplot(614)
            plot(x,Cu,'r',x,Cc,'b','linewidth',2);title('Cu,Cc')
            set(gca,'xlim',[p.xplotmin p.xplotmax])
            
            subplot(615)
            plot(x,Sx,'b','linewidth',2);title('Sx')
            set(gca,'xlim',[p.xplotmin p.xplotmax])
            
            subplot(616)
            plot(x,zb-zb0,'linewidth',2);title('zb-zb0')
            set(gca,'xlim',[p.xplotmin p.xplotmax])
            grid on
            drawnow
            fname=['duna',num2str(round(u10)),'_',num2str(round(phi0*180/pi)),'_',num2str(1000+it),'.jpg']
            print('-djpeg',fname)
            fname=['duna',num2str(round(u10)),'_',num2str(round(phi0*180/pi)),'_',num2str(1000+it)]
            save(fname,'x','u10x','v10x','umagx','u10t0','Cu','Cc','Sx','dzbdt','zb','zbprev','p1','phi','hv','cv')
        end
    end
    
end
p.dt=dt0;