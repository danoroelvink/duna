%% Obter dados de vento com possivel transporte sedimentar
load Vento_Faro2009_2018.txt
load VentoFaro2009_20182.mat %introduzi direcção ventos 2014 com base no site https://www.wunderground.com/history/airport/LPFR/2015/1/2/WeeklyHistory.html?req_city=Faro&req_state=08&req_statename=Portugal&reqdb.zip=00000&reqdb.magic=5&reqdb.wmo=08554
load b0911.txt

direcao_correc= VentoFaro2009_20182(:,3);
Vento_Faro2009_2018dircorrect=[ Vento_Faro2009_2018(:,1:7) direcao_correc Vento_Faro2009_2018(:,9)];
%save('Vento_Faro2009_2017dircorrect.txt','Vento_Faro2009_2017dircorrect','-ascii')

%% selecionar pela direção do vento, neste caso foi entre 127 e 307
ind_th_24=find(Vento_Faro2009_2018dircorrect(:,8)>127 & Vento_Faro2009_2018dircorrect(:,8)<307); 
dir2 = Vento_Faro2009_2018dircorrect(ind_th_24,8);
vento2 = Vento_Faro2009_2018dircorrect(ind_th_24,7);
data2 = Vento_Faro2009_2018dircorrect(ind_th_24,1:6);
year2 = Vento_Faro2009_2018dircorrect(ind_th_24,9);
%% passar vento do aeroporto para a praia com a formula o cftool
ventopraia2009_2018 = vento2*0.979+0.671;

Vento_Faro2009_2018dir = [data2 ventopraia2009_2018 dir2 year2];
%save('Vento_Faro2009_2017dir.txt','Vento_Faro2009_2017dir','-ascii')

%% Compute critical wind velocity as a function of D50 grain size,
% slope angle and moisture content
rho=1.25;
rhos=2650;
g=9.81;
A=0.085;
kappa=0.4;
D=0.000484;  %CASO DAS DUNAS NA PRAIA DE FARO 0.0005
z0=D/30;
Dref = 0.00025;
z=10;
tanphi=tan(35*pi/180);
K =4.2;
phi0=230;
C=1.5;%C=2.8;

ustart=A*sqrt(g*D*((rhos-rho)/rho));

u10t0=ustart/kappa.*log(z./z0); %VELOCIDADE CRITICA
%% selecionar pela critical_velocity
ind_th_2=find(Vento_Faro2009_2018dir(:,7)>u10t0); 
dirvento_final= Vento_Faro2009_2018dir(ind_th_2,8);
vento_final = Vento_Faro2009_2018dir(ind_th_2,7);
data_final = Vento_Faro2009_2018dir(ind_th_2,1:6);
year_final = Vento_Faro2009_2018dir(ind_th_2,9);
%% cos(dir2rad(dirvento-217)

dirvento_final(:,2)=abs(dirvento_final(:,1)-phi0);
dirvento_final(:,3)=deg2rad(dirvento_final(:,2));
dirvento_final(:,4)=cos(dirvento_final(:,3));
dirvento2=dirvento_final(:,4);

%% Compute Fc
Fc=4.38*vento_final(:,7)-8.23;

%% Compute teorical sediment transport - Bagnold 1937 artigo Bas
trans = C*(kappa/(log(z/z0))).^3*600*sqrt(D/Dref)*(rho/g)*(vento_final-ustart).^3;
q=trans/(rhos*0.6).*dirvento2;

%% compare Fc with F (backshore width)
F=b0911;
if Fc <= F
    q=q;
else
    q=q*sin(pi/2*F/Fc);
end

