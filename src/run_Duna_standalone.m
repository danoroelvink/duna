%% Duna: process-based model for aeolian sand transport
%-------------------------------------------------------
% Duna simulates the evolution of dunes as a result of input conditions (e.g., wind, vegetation cover) 
% and of the dynamic interactions and feedbacks between morphology, flow and vegetation. 
% The computational flow in Duna is as follows: 
% a) calculation of the wind velocity field along the beach-dune topographic profile, 
% b) reduction of wind velocity in vegetated areas of the dune due to plant-wind interaction, 
% c) critical shear stress calculation along the profile, accounting for moisture and slope effects for sand transport, 
% d) definition of instantaneous sediment concentrations in the domain, accounting for additional transport limiting parameters (fetch and armoring), and 
% e) update of the dune topography (solution of 1D advection and mass balance) and of the vegetation cover.
%-------------------------------------------------------
%% Model parameterizations are given in: 
% Roelvink and Costas (2019), Coupling nearshore and aeolian processes: XBeach and Duna process-based models. Environmental Modelling & Software, 115, 98–112. https://doi.org/10.1016/J.ENVSOFT.2019.02.010
%% The present version of Duna includes a standalone version that includes parameterisation of the Fetch effect and grid rotation
% Details of the changes in this version can be found in the article 'Exploring controls on coastal dune growth through a simplified model' published in JGR Earth Surface
%% Example application: Faro Beach simulation (experiment Morpho2)
% Prepared for example inputs from a dune profile in Ria Formosa (S. Portugal), in Faro Beach. 
% Profile is near P1 in Costas et al. (2020) https://doi.org/10.1016/j.geomorph.2020.107435 
% Wind data from IPMA, Faro Airport station (https://www.ipma.pt/en/otempo/obs.superficie)
clear all
close all
clc

load('.\inputs\wind_data.mat')
load('.\inputs\dune_profile.mat')

figure
subplot(3,1,1:2)
plot(x.init,z.init,'-k')
hold on
plot(x.final,z.final,'-r')
axis tight
grid on
xlabel('cross-shore distance (m)')
ylabel('elevation (m)')
legend('initial profile','final profile')
subplot(3,1,3)
plot(W.date,W.speed,'-b')
datetick('x','mm/yy')
ylabel('Wind velocity (m/s)')
axis tight
grid on

%% Inputs 
% OffshoreBound='closed';   %no supply considered, not suitable for standalone Duna runs
OffshoreBound='open';       %infinite sand supply from the sea
tide='yes';                 % consider semidiurnal tides with (M2amp, M2freq: amplitude and frequency defined in duna_input)
runname='TestRun';          % name used for the output file
ProfilePrint='yes';         % switch to plot profiles every dtplot days
dtplot=2.0;                 % profile plot step in days
% other inputs defined in duna_input.mat
%%
% ------------ computational grid ------------------
dx=0.5; %fixed horizontal step 
xgrid=dx*ceil(min(x.init)/dx):dx:max(x.init);
zgrid=interp1(x.init,z.init,xgrid,'linear');
phibeach=217;                   % profile direction to N

% _________________________________________________________________
zlee=6.8;%max(zgrid);
xlee=-55;
ztoe=5.5;
zcrest=6.;
% ilee=find(zgrid>=zlee,1,'first');
ilee=find(xgrid>=xlee,1,'first');
icrest=find(zgrid>=zcrest,1,'last');
itoe=find(zgrid>=ztoe,1,'first');
xlee=xgrid(ilee);
xtoe=xgrid(itoe);
xcrest=xgrid(icrest);

ntot=length(W.date)-1;     %simulation steps --> Run all wind timeseries

im=length(xgrid);
hvdata=zeros(im,1);
cvdata=zeros(im,1);
Cc=zeros(im,1);
Cu=zeros(im,1);
CuFc=zeros(im,1);
zb=zgrid;
zbprev = zb;
ztdata=zeros(im,ntot);
hvtdata=zeros(im,ntot);
cvtdata=zeros(im,ntot);
Cudata=zeros(im,ntot);
CuFcdata=zeros(im,ntot);
Ccdata=zeros(im,ntot);
Tau=zeros(im,ntot);

% dats=datenum(wind(:,1:6));
% figure
% plot(dats,wind(:,7),'.r')

durats=diff(W.date)*86400.;


istorm=zeros(ntot,1);

Emma1=datenum('28/02/2018 04:00','dd/mm/yyyy HH:MM');
Emma2=datenum('05/03/2018 20:00','dd/mm/yyyy HH:MM');
iEmm=find(W.date>=Emma1&W.date<=Emma2);
istorm(iEmm,:)=1;


iplot=zeros(ntot);
time=W.date-W.date(1);

switch ProfilePrint
    case 'yes'
ipr=find(rem(time,dtplot)==0);
iplot(ipr)=1;
end


%%

for n=1:ntot
zbb=zb;
durat=durats(n);
dir=deg2rad(W.dir(n)-phibeach);
u10=W.speed(n);
durat_nowind=0.0;

if n==iEmm(end)+1
%changing profile    
ijk=find(zb>=5.2,1,'first');
xx1=[x.init(1),xgrid(ijk)];
zz1=[z.init(1),zb(ijk)];
zbb(2:ijk-1)=interp1(xx1,zz1,xgrid(2:ijk-1));
end
zgrid=zbb;
xrunup=min(xgrid);%0.0;

hv = hvdata;
cv = cvdata;


storm=istorm(n,1);    
[zb,hv,cv,Cc,Cu,CuFc,sandtrans]= duna (xgrid,zbb,itoe,icrest,ilee,durat,dir,u10,zbprev,hv,cv,xrunup,durat_nowind,OffshoreBound,tide,storm);

hvdata=hv;
cvdata=cv;

dzb=(zb-zbb)';
zbprev=zbb; %old timestep
zbb=zb;



ztdata(:,n)=zb;
hvtdata(:,n)=hv;
cvtdata(:,n)=cv;
Cudata(:,n)=Cu;
CuFcdata(:,n)=CuFc;
Ccdata(:,n)=Cc;
switch sandtrans
    case 'yes'
load ('kroy.mat')
Tau(:,n)=tau_over_tau0;
    case 'no'
Tau(:,n)=zeros(1,im);    
end

% updating the locations of the toe and lee after morpho changes
ilee=find(zgrid>=zlee,1,'first');
itoe=find(zgrid>=ztoe,1,'first');
xlee=xgrid(ilee);
xtoe=xgrid(itoe);

switch num2str(iplot(n))
    case '1'
figure(111)
plot(x.init,z.init,'-k')
hold on
plot(xgrid,zb,'-r')
grid on 
axis tight 
legend('initial',['step=',num2str(n)],'location','southeast')
text(0.98*xgrid(1), 0.98*max(max(ztdata)),{'simulated profile';'at day:';[num2str(time(n)), ' from ', num2str(time(ntot),'%.0f')]},'VerticalAlignment', 'top','Backgroundcolor','k','color','w','FontWeight','bold')
% ***** plotting simulated vegetation *****
ii=find(cv>0);
for k=1:length(ii)
    j=ii(k);
plot([xgrid(j) xgrid(j)],[zb(j) zb(j)+hv(j)],'-','color',[153 204 51]/255)
end
% ******************************************
hold off
end


end

    save([runname,'.mat'],'xgrid', 'zgrid', 'hvtdata', 'cvtdata', 'ztdata', 'Cudata', 'CuFcdata', 'Ccdata', 'Tau')

%% calculating simulated and obvserved volumes and RMSE 
zz1=5.5; %seaward-most elevation to calculate dune volume
zz2=4.7; %landward-most elevation to calculate dune volume

zzo2=ztdata(:,end);
zzo1=ztdata(:,1);
io1=find(zzo1>=5.5,1,'first');
dvol=trapz(xgrid(ioo:im),zzo2(ioo:im))-trapz(xgrid(ioo:im),zzo1(ioo:im));


j1a=find(z.init>=zz1,1,'first');
j1b=find(z.init(j1a:end)>=zz2,1,'last');
j1b=j1b+j1a;
j2a=find(z.final>=zz1,1,'first');
j2b=find(z.final(j2a:end)>=zz2,1,'last');
j2b=j2b+j2a;
dvolo=trapz(x.final(j2a:j2b),z.final(j2a:j2b))-trapz(x.init(j1a:j1b),z.init(j1a:j1b)); 

% Interpolating measured profile on the grid
zmeas=interp1(x.final,z.final,xgrid);
Err=zmeas(ioo:im)-zzo2(ioo:im)';
Err=Err.^2.0;
Nn=im-ioo+1;
RMSE=sqrt(sum(Err)/Nn);

resVol={strcat('dV sim=',num2str(dvol,3));strcat('dV obs=',num2str(dvolo,3));strcat('dV sim-obs=',num2str(dvol-dvolo,3));strcat('RMSE=',num2str(RMSE,3),'m')};
disp('---------- End of Duna Simulation ------------')
disp(['dV sim=',num2str(dvol)])
disp(['dV obs=',num2str(dvolo)])
disp(['dV sim-obs=',num2str(dvol-dvolo)])
disp(['RMSE=',num2str(RMSE),' m'])


figure
lin1=plot(xgrid,ztdata(:,1),'-k');
hold on
lin2=plot(xgrid,ztdata(:,n),'-r');
plot(xgrid,ztdata(:,iEmm(1)-1),'-b');
plot(xgrid,ztdata(:,iEmm(end)+1),'-g');
ii=find(cvtdata(:,n)>0);
grr=[153 204 51]/255;
for k=1:length(ii)
    j=ii(k);
plot([xgrid(j) xgrid(j)],[ztdata(j,n) ztdata(j,n)+hvtdata(j,n)],'-','color',grr,'linewidth',cvtdata(j,n))
end
lin3=plot(xgrid,zmeas,'--c','linewidth',2);
grid on
axis tight
hold off
ylabel('elevation [m]')
xlabel('cross-shore distance [m]')
text(0.98*xgrid(1), max(max(ztdata)),resVol,'VerticalAlignment', 'top','Backgroundcolor','k','color','w','FontWeight','bold')
legend([lin1 lin2 lin3],'initial zb','final zb','measured','location','northeast')   
