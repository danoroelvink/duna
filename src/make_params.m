function []=make_params(modeldir,nx,event,no,Hm0,Tp,dir,zs0,lsgrad,tstop,morfac)
switch event
    case 'tsunami'
        swave=0;
        nonh=1;
        instat='ts_nonh';
        front='nonh_1d';
        morfac=1;
        morfacopt=0;
        tintg=10;
        morstart=0;
        breakform=2;
        gamma=0.75;
        taper=60;
        gammax=.8;
        wavint=60;
    case 'moderate'
        swave=1;
        nonh=0;
        instat='stat_table';
        front='abs_1d';
        morfacopt=0;
        tintg=60;
        morstart=600;
        tstop=tstop+morstart;
        breakform=2;
        gamma=0.75;
        taper=60;
        gammax=.8;
        wavint=60;
        bermslope=0.15;
    case 'storm'
        swave=1;
        nonh=0;
        instat='jons_table';
        front='abs_1d';
        morfacopt=0;
        tintg=5;
        morstart=600;
        tstop=tstop+morstart;    
        breakform=3;
        gamma=0.55;
        taper=60;
        gammax=2;
        wavint=0;
        bermslope=0.15;
    case 'moderate_nonh'
        swave=0;
        nonh=1;
        instat='jons_table';
        front='nonh_1d';
        morfacopt=0;
        lwt= 1;
        betad= 1;
        tintg=5;
        morstart=600;
        tstop=tstop+morstart;    
        taper=60;
        wavint=0;
end
        
fid=fopen([modeldir,'/','params.txt'],'w');
fprintf(fid, '%s', ['ships = 0']);
fprintf(fid, '\n%s', ['swave = ',num2str(swave)]);
fprintf(fid, '\n%s', ['lwave = 1']);
fprintf(fid, '\n%s', ['sedtrans = 1']);
fprintf(fid, '\n%s', ['morphology = 1']);
fprintf(fid, '\n%s', ['avalanching = 1']);
fprintf(fid, '\n%s', ['gwflow = 1']);
fprintf(fid, '\n%s', ['gwnonh = 1']);
fprintf(fid, '\n%s', ['kx = 0.01']);
fprintf(fid, '\n%s', ['ky = 0.01']);
fprintf(fid, '\n%s', ['kz = 0.01']);
fprintf(fid, '\n%s', ['tsfac= 0.2'])
fprintf(fid, '\n%s', ['secorder = 0']);
fprintf(fid, '\n%s', ['nonh    = ',num2str(nonh)]);
fprintf(fid, '\n%s', ['bermslope    = ',num2str(bermslope)]);
fprintf(fid, '\n%s', [' ']);


fprintf(fid, '\n%s', ['nx=',num2str(nx)]);
fprintf(fid, '\n%s', ['ny=0']);
fprintf(fid, '\n%s', ['vardx=1']);
fprintf(fid, '\n%s', ['depfile=','profile.dep']);
fprintf(fid, '\n%s', ['xfile=','x.dep']);
fprintf(fid, '\n%s', ['posdwn=0 ']);
fprintf(fid, '\n%s', ['thetanaut=0']);
fprintf(fid, '\n%s', ['thetamin=-90']);
fprintf(fid, '\n%s', ['thetamax= 90']);
fprintf(fid, '\n%s', ['dtheta=180']);
fprintf(fid, '\n%s', [' ']);

if ~strcmp(event,'moderate_nonh')
fprintf(fid, '\n%s', ['break=',num2str(breakform)]);
fprintf(fid, '\n%s', ['gamma=',num2str(gamma)]);
fprintf(fid, '\n%s', ['gammax=',num2str(gammax)]);
fprintf(fid, '\n%s', ['wavint=',num2str(wavint)]);
fprintf(fid, '\n%s', ['n=10']);
fprintf(fid, '\n%s', ['fw=0.01']);
fprintf(fid, '\n%s', ['lwt=1']);
fprintf(fid, '\n%s', ['betad=1']);
fprintf(fid, '\n%s', ['turbadv=lagrangian']);
fprintf(fid, '\n%s', [' ']);
end

fprintf(fid, '\n%s', ['front = ',front]);
fprintf(fid, '\n%s', ['left  = neumann']);
fprintf(fid, '\n%s', ['right = neumann']);
fprintf(fid, '\n%s', ['back  = abs_1d']);
fprintf(fid, '\n%s', [' ']);

if strcmp(event,'moderate')
    fprintf(fid, '\n%s', ['tideloc = 1']);
    fprintf(fid, '\n%s', ['zs0file = tide.txt']);
    fprintf(fid, '\n%s', ['facsl= 0.5'])
end
if strcmp(event,'moderate_nonh')
    fprintf(fid, '\n%s', ['tideloc = 1']);
    fprintf(fid, '\n%s', ['zs0file = tide.txt']);
end
if strcmp(event,'storm')
    fprintf(fid, '\n%s', ['tideloc = 1']);
    fprintf(fid, '\n%s', ['zs0file = tide.txt']);
end
fprintf(fid, '\n%s', ['CFL = 0.3']);
fprintf(fid, '\n%s', ['taper =',num2str(taper)]);
fprintf(fid, '\n%s', ['tstop = ',num2str(tstop)]);
fprintf(fid, '\n%s', ['tintg = ',num2str(tintg)]);
fprintf(fid, '\n%s', ['morfac= ',num2str(morfac)]);
fprintf(fid, '\n%s', ['morfacopt= ',num2str(morfacopt)]);
fprintf(fid, '\n%s', ['morstart= ',num2str(morstart)]);
fprintf(fid, '\n%s', ['lsgrad= ',num2str(lsgrad)]);
fprintf(fid, '\n%s', ['facSk= 0.15']);
fprintf(fid, '\n%s', ['facAs= 0.35']);
fprintf(fid, '\n%s', ['tsfac= 0.2'])
fprintf(fid, '\n%s', ['D50= 0.0005']);
fprintf(fid, '\n%s', ['D90= 0.0008']);

fprintf(fid, '\n%s', ['waveform= ruessink_vanrijn']);
fprintf(fid, '\n%s', ['turb = wave_averaged']);
fprintf(fid, '\n%s', ['bedfriction = manning']);
fprintf(fid, '\n%s', ['bedfriccoef = 0.025']);
fprintf(fid, '\n%s', ['maxbrsteep = 0.3']);
fprintf(fid, '\n%s', [' ']);

fprintf(fid, '\n%s', ['depthscale=1']);
fprintf(fid, '\n%s', [' ']);

fprintf(fid, '\n%s', ['zs0          = ',num2str(zs0)]);
fprintf(fid, '\n%s', [' ']);

fprintf(fid, '\n%s', ['instat       = ',instat]);
fprintf(fid, '\n%s', ['bcfile       = jonswap.txt']);
fprintf(fid, '\n%s', ['random       = 0']);
fprintf(fid, '\n%s', [' ']);

% fprintf(fid, '\n%s', ['nrugauge     =   1']);
% fprintf(fid, '\n%s', ['-5 5']);
% fprintf(fid, '\n%s', [' ']);

fprintf(fid, '\n%s', ['outputformat = netcdf']);
fprintf(fid, '\n%s', ['ncfilename=xboutput.nc']);
fprintf(fid, '\n%s', [' ']);

fprintf(fid, '\n%s', ['nglobalvar = 8']);
fprintf(fid, '\n%s', ['zb']);
fprintf(fid, '\n%s', ['zs']);
fprintf(fid, '\n%s', ['H']);
fprintf(fid, '\n%s', ['ue']);
fprintf(fid, '\n%s', ['ve']);
fprintf(fid, '\n%s', ['Susg']);
fprintf(fid, '\n%s', ['Svsg']);
fprintf(fid, '\n%s', ['sedero']);
fclose(fid);

fid=fopen([modeldir,'/','jonswap.txt'],'w');
for i=1:ceil((tstop+1)/3600)
    fprintf(fid, ' %s \n', [num2str(Hm0),' ',num2str(Tp),' ',num2str(dir), ...
        ' 3.3 1000 ',num2str(3600),' 1']);
end
fclose(fid);
end
