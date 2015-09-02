if strcmp(lake,'Crater')==1 || strcmp(lake,'Superior')==1 || strcmp(lake,'Erie')==1 || strcmp(lake,'Michigan')==1 || strcmp(lake,'Huron')==1 || strcmp(lake,'Ontario')==1  || strcmp(lake,'Baikal')==1 || strcmp(lake,'Sparkling')==1 || strcmp(lake,'Crystal')==1 
    min_maxrad=270; max_maxrad=320; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=0;   max_minrad=50;  delta_minrad=max_minrad-min_minrad;
elseif strcmp(lake,'Tahoe')==1
    min_maxrad=325; max_maxrad=410; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=30;   max_minrad=110;  delta_minrad=max_minrad-min_minrad;
elseif strcmp(lake,'Garda')==1
    min_maxrad=270; max_maxrad=330; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=0;   max_minrad=75;  delta_minrad=max_minrad-min_minrad;
elseif strcmp(lake,'Mara')==1
    min_maxrad=250; max_maxrad=300; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=0;   max_minrad=50;  delta_minrad=max_minrad-min_minrad;
elseif strcmp(lake,'Neusiedl')==1  % see http://www.zamg.ac.at/fix/klima/oe71-00/klima2000/klimadaten_oesterreich_1971_frame1.htm
    min_maxrad=220; max_maxrad=270; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=10;   max_minrad=50;  delta_minrad=max_minrad-min_minrad;
elseif strcmp(lake,'Constance')==1  % see http://www.zamg.ac.at/fix/klima/oe71-00/klima2000/klimadaten_oesterreich_1971_frame1.htm
% and also http://www.meteoswiss.admin.ch/web/it/clima/clima_della_svizzera/norma_1961_90.html
    min_maxrad=180; max_maxrad=220; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=10;   max_minrad=50;  delta_minrad=max_minrad-min_minrad;
elseif strcmp(lake,'Zurich')==1 
    min_maxrad=200; max_maxrad=280; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=10;   max_minrad=50;  delta_minrad=max_minrad-min_minrad;
elseif strcmp(lake,'Biel')==1
    min_maxrad=200; max_maxrad=280; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=10;   max_minrad=50;  delta_minrad=max_minrad-min_minrad;
elseif strcmp(lake,'Balaton')==1 % http://w3.georgikon.hu/tanszekek/Meteor/docs/aa9.pdf and also http://2ndwinterlimnology.igb-berlin.de/Presentations/Voros.pdf
    min_maxrad=310; max_maxrad=360; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=25;   max_minrad=75;  delta_minrad=max_minrad-min_minrad;
elseif strcmp(lake,'SFrancisco')==1 % A Cubic Mile of Oil: Realities and Options for Averting the Looming Global ... pag 197-198
    min_maxrad=280; max_maxrad=340; delta_maxrad=max_maxrad-min_maxrad;
    min_minrad=30;   max_minrad=100;  delta_minrad=max_minrad-min_minrad;
end

maxrad=min_maxrad:delta_maxrad/n:max_maxrad;
minrad=min_minrad:delta_minrad/n:max_minrad;

min_rs=0.04; max_rs=0.20; delta=max_rs-min_rs;  % shortwave reflectivity (albedo)
rs=min_rs:delta/n:max_rs;

for i=1:n+1
    for j=1:n+1
        Rad_a(i,j)=(maxrad(i)-minrad(j))/2;
        Rad_b(i,j)=minrad(j)+Rad_a(i);
    end
end

minRad_a=min(min(Rad_a));    maxRad_a=max(max(Rad_a));
minRad_b=min(min(Rad_b));    maxRad_b=max(max(Rad_b));

min_Temperatura_rif=0; max_Temperatura_rif=30; delta=max_Temperatura_rif-min_Temperatura_rif;
Temperatura_rif=min_Temperatura_rif:delta/n:max_Temperatura_rif;

min_TW=0; max_TW=30; delta=max_TW-min_TW;
TW=min_TW:delta/n:max_TW;

min_aa=0.6; max_aa=0.9; delta=max_aa-min_aa;
aa=0.97*(min_aa:delta/n:max_aa);

bb=0.97;

min_alphac=3; max_alphac=15; delta=max_alphac-min_alphac;
alphac=min_alphac:delta/n:max_alphac;

Kelvin0=273.15 + Temperatura_rif;	

min_Delta_alphac=0.1; max_Delta_alphac=15; delta=max_Delta_alphac-min_Delta_alphac;
Delta_alphac=min_Delta_alphac:delta/n:max_Delta_alphac;

min_ea=5; max_ea=15; delta=max_ea-min_ea;
ea=min_ea:delta/n:max_ea;

min_Delta_ea=0.1; max_Delta_ea=10; delta=max_Delta_ea-min_Delta_ea;
Delta_ea=min_Delta_ea:delta/n:max_Delta_ea;

min_DT=0; max_DT=30; delta=max_DT-min_DT;
DT=min_DT:delta/n:max_DT;

for i=1:n+1
	for j=1:n+1
		ew_par(i,j)=6.112*exp(Temperatura_rif(i)*17.67/(Temperatura_rif(i)+243.5)); %* ...	
		%(1+17.67*243.5/(Temperatura_rif(i)+243.5)^2*(TW(j)-Temperatura_rif(i)));
	end
end
ew=[max(0,min(min(ew_par))) max(max(ew_par))];

alphae=alphac/0.61;          

Delta_alphae = Delta_alphac;

% questi sono noti, non serve modificarli
rho=1000;
cp=4186;
s_Boltzman=5.67E-8;

