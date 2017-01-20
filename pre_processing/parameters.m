min_maxrad=200; max_maxrad=450;  delta_maxrad=max_maxrad-min_maxrad;
min_minrad=0;   max_minrad=250;  delta_minrad=max_minrad-min_minrad;

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

% constants
rho=1000;
cp=4186;
s_Boltzman=5.67E-8;

