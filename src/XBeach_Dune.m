%% XBeach_D
%  Program to simulate combined effects of moderate waves, storms,
%  aeolian transport and tsunamis on profile evolution and stratigraphy

%% General input

xb = XB_dune_input()
if ~isfield(xb,'outputopt')
    xb.outputopt=1;
end
if ~isfield(xb,'windfactor')
    xb.windfactor=1;
end
if xb.restart==0
    clear all; close all
    xb = XB_dune_input()
    if ~isfield(xb,'outputopt')
        xb.outputopt=1;
    end
    if ~isfield(xb,'windfactor')
        xb.windfactor=1;
    end
   addpath (xb.functions_path);
    %% Read initial profile
    profname=xb.profname;
    xz=load(profname);
    dxmother=xb.dxmother;
    xdata1=-xz(end:-1:1,1);
    zdata1=xz(end:-1:1,2);
    xdata=[xdata1(1):dxmother:xdata1(end)];
    zdata=interp1(xdata1,zdata1,xdata);
    ztdata(1,:)=zdata;
    hvdata=zeros(size(xdata));
    hvtdata(1,:)=hvdata;
    cvdata=zeros(size(xdata));
    cvtdata(1,:)=cvdata;
    Sydata=zeros(size(xdata));
    Sytdata(1,:)=Sydata;
    %   zbprevdata=zeros(size(xdata));
    zbprevdata=zdata;
    fname=[num2str(1000),'.mat'];
    if xb.outputopt==2
        save(fname,'xdata','zdata','hvdata','cvdata','Sydata');
    end
else
    try
        fname=xb.restartfile;
    catch
        fname=[num2str(1000+xb.restart),'.mat'];
    end
    load(fname);
    zdata=ztdata(xb.restart,:);
    hvdata=hvtdata(xb.restart,:);
    cvdata=cvtdata(xb.restart,:);
    zbprevdata=zdata;
end

%% Per event type read event conditions table
cond_tsunami    = load(xb.cond_tsunami);
cond_moderate   = load(xb.cond_moderate);
cond_storm      = load(xb.cond_storm);
cond_dune       = load(xb.cond_dune);

%% Read list of events to simulate sequentially
%  very simple coding:
%  event_type event_number
%  e.g.
%  tsunami 1
%  moderate 3
%  dune 1
%  moderate 5
%  dune 1
%  storm 5
fid=fopen(xb.events);
i=1
first_dune=1;
t(1)=0;
xrunup=-1e10;
dur_nowind=0;
while 1
    %% Read next event
    line = fgetl(fid);
    if ~ischar(line), break, end
    [str1 str2]=strtok(line);
    event{i}=str1;
    no(i)=str2num(str2);
    switch event{i}
        case 'moderate'
            u10in(i)   = nan;
            Hm0   =cond_moderate(no(i),1);
            Hm0t(i)=Hm0;
            tstop =cond_moderate(no(i),4);
            morfac=cond_moderate(no(i),5);
            dtmor=tstop*morfac;
            dur_nowind=dur_nowind+dtmor;
            durat_nowind(i)=nan;
            t(i+1)=t(i)+tstop*morfac;
            twd(i+1)=nan;
         case 'storm'
            u10in(i)   = nan;
            tstop =cond_storm(no(i),4);
            morfac=cond_storm(no(i),5);
            Hm0   =cond_storm(no(i),1);
            Hm0t(i)=Hm0;
            dtmor=tstop*morfac;
            dur_nowind=dur_nowind+dtmor;
            durat_nowind(i)=nan;
            t(i+1)=t(i)+tstop*morfac;
            twd(i+1)=nan;
            
        case 'moderate_nonh'
            u10in(i)   = nan;
            Hm0   =cond_moderate_nonh(no(i),1);
            Hm0t(i)=Hm0;
            tstop =cond_moderate_nonh(no(i),4);
            morfac=cond_moderate_nonh(no(i),5);
            dtmor=tstop*morfac;
            dur_nowind=dur_nowind+dtmor;
            durat_nowind(i)=nan;
            t(i+1)=t(i)+tstop*morfac;
            twd(i+1)=nan;
            
        case 'dune'
            u10 =cond_dune(no(i),1)*xb.windfactor;
            phi0=cond_dune(no(i),2);
            durat =cond_dune(no(i),3);
            Hm0t(i) = nan;
            durat_nowind(i)=max(dur_nowind-durat,0);
            dur_nowind=0;
            u10in(i)=u10;
            t(i+1)=t(i);%+durat;
            twd(i+1)=t(i)+durat;
            phi0=deg2rad(abs(phi0-xb.phibeach));
            
    end
    i=i+1
end
for i=xb.restart+1:length(Hm0t)
    if i>0%size(ztdata,1)-1
        %% Prepare profile for simulation
        for j=1:length(xb.eventname)
            if strcmp(event{i},xb.eventname{j})
                [xgrid zgrid] = makeprofile(xdata,zdata, ...
                    xb.zbmin(j),xb.zbmax(j),xb.dxmin(j),xb.dxmax(j),xb.phibeach);
                dzb=zeros(size(xgrid));
                Sy=zeros(size(xgrid));
            elseif strcmp(event{i},'stop')
                break
            end
        end
        %% Prepare input for simulation
        modeldir='.';
        nx=length(xgrid)-1;
        switch event{i}
            case 'tsunami'
                if xb.xbeach == 1
                    H = cond_tsunami(no(i),1);
                    zs0=cond_tsunami(no(i),2);
                    [t_ts,eta_ts]=solitary (H,100);
                    tstop = t_ts(end);
                    make_params(modeldir,nx,xb.phibeach,event{i},no(i),0,0,0,zs0,0,tstop,1)
                end
            case 'moderate'
                if xb.xbeach == 1
                    Hm0   =cond_moderate(no(i),1);
                    Tp    =cond_moderate(no(i),2);
                    dir   =cond_moderate(no(i),3);
                    tstop =cond_moderate(no(i),4);
                    morfac=cond_moderate(no(i),5);
                    zs0   = maketide2(xb.MLWS,xb.MHWS,tstop);
                    make_params(modeldir,nx,xb.phibeach,event{i},no(i),Hm0,Tp,dir,zs0,xb.lsgrad,tstop,morfac)
                    dtmor=tstop*morfac;
                end
                
            case 'storm'
                if xb.xbeach == 1
                    Hm0   =cond_storm(no(i),1);
                    Tp    =cond_storm(no(i),2);
                    dir   =cond_storm(no(i),3);
                    tstop =cond_storm(no(i),4);
                    morfac=cond_storm(no(i),5);
                    % zs0   =cond_storm(no(i),6);
                    zs0   = maketide2(xb.MLWS,xb.MHWS,tstop);
                    make_params(modeldir,nx,xb.phibeach,event{i},no(i),Hm0,Tp,dir,zs0,xb.lsgrad,tstop,morfac)
                    dtmor=tstop*morfac;
                end
                
            case 'moderate_nonh'
                if xb.xbeach == 1
                    Hm0   =cond_moderate_nonh(no(i),1);
                    Tp    =cond_moderate_nonh(no(i),2);
                    dir   =cond_moderate_nonh(no(i),3);
                    tstop =cond_moderate_nonh(no(i),4);
                    morfac=cond_moderate_nonh(no(i),5);
                    zs0   = 0;
                    lsgrad=0;
                    make_params(modeldir,nx,xb.phibeach,event{i},no(i),Hm0,Tp,dir,zs0,lsgrad,tstop,morfac)
                    dtmor=tstop*morfac;
                end
            case 'dune'
                u10 =cond_dune(no(i),1)*xb.windfactor;
                durat =cond_dune(no(i),3);
                dir=cond_dune(no(i),2);
                dir=deg2rad(abs(dir-217));
        end
        
        %% Run simulation
        switch event{i}
            case 'tsunami'
                if xb.xbeach==1
                   !run
                end
            case 'moderate'
                if xb.xbeach==1
                   !run
                end
            case 'storm'
                if xb.xbeach==1
                   !run
                end
            case 'moderate_nonh'
                if xb.xbeach==1
                   !run
                end
            case 'dune'
                if xb.duna==1
                    hv = interp1(xdata,hvdata,xgrid);
                    cv = interp1(xdata,cvdata,xgrid);
                    zbprev = interp1(xdata,zbprevdata,xgrid);
                    [zb,hv,cv]= duna (xgrid,zgrid,durat,dir,u10,zbprev,hv,cv,xrunup,durat_nowind(i));
                    dzb=(zb-zgrid)';
                    zbprev=zb;
                end
        end
        if 1
            %% Extract end profile and transports
            switch event{i}
                case 'tsunami'
                    if xb.xbeach==1
                        %x=nc_varget('xboutput.nc','globalx');
                        %Note: this goes wrong if grid is turned so use xgrid!
                        x=xgrid';
                        zb=nc_varget('xboutput.nc','zb');
                        zblast=squeeze(zb(end,1,:));
                        zbfirst=squeeze(zb(1,1,:));
                        dzb=zblast-zbfirst;
                        Sy=zeros(size(x));
                        xrunup=max(xrunup,getrunup('xboutput.nc'));
                    end
                case 'moderate'
                    if xb.xbeach==1
                        %x=nc_varget('xboutput.nc','globalx');
                        %Note: this goes wrong if grid is turned so use xgrid!
                        x=xgrid';
                        zb=nc_varget('xboutput.nc','zb');
                        zblast=squeeze(zb(end,1,:));
                        zbfirst=squeeze(zb(1,1,:));
                        dzb=zblast-zbfirst;
                        Sy=mean(squeeze(nc_varget('xboutput.nc','Svsg')),1);
                        xrunup=max(xrunup,getrunup('xboutput.nc'));
                        %copyfile('xboutput.nc',['xb',num2str(1000+i),'.nc']);
                    end
                case 'storm'
                    if xb.xbeach==1
                        %x=nc_varget('xboutput.nc','globalx');
                        %Note: this goes wrong if grid is turned so use xgrid!
                        x=xgrid';
                        zb=nc_varget('xboutput.nc','zb');
                        zblast=squeeze(zb(end,1,:));
                        zbfirst=squeeze(zb(1,1,:));
                        dzb=zblast-zbfirst;
                        Sy=mean(squeeze(nc_varget('xboutput.nc','Svsg')),1);
                        xrunup=max(xrunup,getrunup('xboutput.nc'));
                        %copyfile('xboutput.nc',['xb',num2str(1000+i),'.nc']);
                    end
                case 'dune'
                    if xb.duna==1
                        x=xgrid';
                        Sy=zeros(size(x));
                        %              zblast=zb';
                        inmodel=xdata>=xgrid(1)&xdata<=xgrid(end);
                        hvdata(inmodel)=interp1(x,hv,xdata(inmodel),'linear','extrap');
                        cvdata(inmodel)=interp1(x,cv,xdata(inmodel),'linear','extrap');
                        zbprevdata(inmodel)=interp1(x,zbprev,xdata(inmodel),'linear','extrap');
                        xrunup=-1e10;
                    end
            end
            %% Interpolate back to mother profile
            if max(abs(dzb))>5
                disp ('maximum bed level change > 5 m')
                return
            end
            inmodel=xdata>=xgrid(3)&xdata<=xgrid(end);
            zdata(inmodel)=zdata(inmodel)+interp1(xgrid,dzb,xdata(inmodel),'linear','extrap');
            ztdata(i+1,:)=zdata;
            hvtdata(i+1,:)=hvdata;
            cvtdata(i+1,:)=cvdata;
            Sydata=zeros(size(xdata));
            Sydata(inmodel)=interp1(xgrid,Sy,xdata(inmodel),'linear','extrap');
            Sytdata(i,:)=Sydata;
            %close all
            if xb.outputopt==1
                %% Update stratigraphy
                tstrat=[0:i];
                [ xtdata,zstrat,time ] = stratigraphy( xdata,ztdata,tstrat );
            end
            if mod(i,xb.plotinterval)==0
                figure(1)
                Hm0=Hm0t(i);
                g=subplot(4,2,[1 3]);
                hold off
                f1=plot(xdata,ztdata,'Color', [0.4,0.4,0.4],'linewidth',0.5);
                hold on
                f2=plot(xdata,ztdata(1,:),'b','linewidth',2);
                f3=plot(xdata,ztdata(end,:),'r','linewidth',2);
                if isnan(u10in(i))
                    annotation('textbox',[0.15 0.715 0.2 0.1],'String',['Hm0 =' num2str(Hm0)],...
                        'LineStyle','none', 'BackgroundColor','w','FontSize',9,'Color','r')
                else
                    annotation('textbox',[0.15 0.715 0.2 0.1],'String',['U =' num2str(u10in(i),'%2.1f')],...
                        'LineStyle','none','BackgroundColor','w','FontSize',9,'Color','r')
                end
                legend([f2 f3],{'Initial','End'},'Location','northwest');
                legend('boxoff');
                set(gca,'Xlim',[-450 0]);
                set(gca,'Ylim',[-5 8]);
                xlabel('Cross-shore distance (m)','FontSize',9);
                ylabel('Elevation (m, MSL)','FontSize',9);
                
                twd1=t(2:end);
                twd2=twd(2:end);
                fd=u10in;
                fd0=zeros(size(fd));
                wind=[1 0.8 0];
                figure(1)
                subplot('position',[0.12, 0.27, 0.8, 0.16]);
                %subplot(4,2,[5 6])
                plot([twd1;twd2],[u10in;u10in],'.-b','linewidth',2);
                fill([twd1;twd2;twd2;twd1;twd1],[fd;fd;fd0;fd0;fd],wind,'linewidth',1);
                if ~isnan(fd(i))
                    hold on;
                    plot([twd2(i),twd2(i)],[0,20],'r','linewidth',2)
                    hold off;
                end
                xmax=max(twd);
                xlim([0 xmax]);
                ylabel({'Wind speed';'(m/s)'},'FontSize',9);
                
                t1=t(1:end-1);
                t2=t(2:end);
                f=Hm0t;
                f0=zeros(size(f));
                figure(1)
                subplot('position',[0.12, 0.1, 0.8, 0.16]);
                wave=[0 0.5 1];
                plot([t1;t2]/24/3600/30.5,[Hm0t;Hm0t],'.-k','linewidth',2);
                fill([t1;t2;t2;t1;t1]/24/3600/30.5,[f;f;f0;f0;f],wave,'linewidth',1)
                if ~isnan(Hm0t(i))
                    hold on;
                    plot([t2(i),t2(i)]/24/3600/30.5,[0,10],'r','linewidth',2)
                    hold off;
                end
                xlim([0 xmax/24/3600/30.5]);
                xlabel('Morfological time (months)');
                ylabel({'Wave height';'(m)'},'FontSize',9);
                
                figure(1)
                k=subplot(4,2,[2 4]);
                pcolor(xtdata,zstrat,time);
                hold on
                plot(xtdata(1,:),zstrat(1,:),'k',xtdata(end,:),zstrat(end,:),'b')
                hold off
                shading flat;
                colormap jet;
                %c=colorbar;
                %cred = 0.05*c;
                xlabel ('x (m)')
                ylabel ('bed level (m)')
                title ('Age of sediment (#)')
                set(gca,'xlim',[-450 0]);
                set(gca,'ylim',[-5 8]);
                fname=[num2str(2000+i),'.jpg'];
                print('-djpeg',fname)
            end
            fname=[num2str(1000+i),'.mat'];
            if xb.outputopt==1
                save(fname,'xdata','xtdata','ztdata','zstrat','hvtdata','cvtdata','Sytdata','time','t','event','no');
                fnamedel=[num2str(1000+i-1),'.mat'];
                try
                    delete(fnamedel);
                end
            else
                save(fname,'xdata','zdata','hvdata','cvdata','Sydata');
            end
        end
    end
end
fclose(fid);










