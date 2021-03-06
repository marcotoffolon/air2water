%       .__       ________                  __
% _____  |__|______\_____  \__  _  _______ _/  |_  ___________
% \__  \ |  \_  __ \/  ____/\ \/ \/ /\__  \\   __\/ __ \_  __ \
%  / __ \|  ||  | \/       \ \     /  / __ \|  | \  ___/|  | \/
% (____  /__||__|  \_______ \ \/\_/  (____  /__|  \___  >__|
%      \/                  \/             \/          \/
% A model to predict Lake Surface Temperature (LST) using air temperature.
% Version 2.0.0 - January 2017
%
% Provided by Sebastiano Piccolroaz and Marco Toffolon
% Department of Civil, Environmental, and Mechanical Engineering, University of Trento (Italy)
% email contacts: s.piccolroaz@unitn.it, marco.toffolon@unitn.it

% PRE-PROCESSING
% Aim: To estimate the a priori range of variation of model parameters

% close all
clear all
clc

disp('This script evaluates the a-priori range of variation of model parameters for temperate lakes')

depth=input('What is the mean depth of the lake for which you want to calculate the a priori range of variation of model parameters?')


set(0,'defaulttextfontname','Times New Roman','defaultaxesfontname','Times New Roman',...
    'defaulttextfontsize',8,'defaultaxesfontsize',8);

n=2;

% call parameters.m
parameters

meanD=[1:9 10:10:100 200:100:1000];    % Mean depth array used to draw the figure
meanD=[meanD depth];
meanD=sort(meanD);

for prof=1:length(meanD)
    min_D = 1+meanD(prof)/1000*50;
    max_D = max(10,meanD(prof));
    
    % Definition: the reactive volume is the volume participating to the heat
    % exchange with the atmosphere. The reactive volume is maximum when the lake
    % is not stratified.
    
    delta=max_D-min_D;
    D=min_D:delta/n:max_D;
    
    pp=rho*cp*D/86400;
    
    % For a detailed derivation of parameters, see Supplementary Material in Piccolroaz (2016).
    
    for i=1:n+1
        for j=1:n+1
            for k=1:n+1
                for l=1:n+1
                    for m=1:n+1
                        for o=1:n+1
                            for p=1:n+1
                                for q=1:n+1
                                    temp=6.112*exp(Temperatura_rif(m)*17.67/ ...
                                        (Temperatura_rif(m)+243.5))*( 1 - 17.67*243.5/ ...
                                        (Temperatura_rif(m)+243.5)^2*Temperatura_rif(m) );
                                    
                                    p1_par(i,j,k,l,m,o,p,q)=( ...
                                        (1-rs(p))*Rad_b(i,i) + ... 
                                        s_Boltzman*Kelvin0(m)^3* ...
                                        (aa(j)-bb)*(273.15-3*Temperatura_rif(m))...
                                        + alphae(k)*( ea(l) - temp ) ...  % note that the sign of this term is wrong in Piccolroaz et al., 2013
                                        )/pp(o);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    p1=[min(min(min(min(min(min(min(min(p1_par)))))))) max(max(max(max(max(max(max(max(p1_par))))))))];
    p2=[(4*s_Boltzman*aa(1)*Kelvin0(1)^3 + alphac(1))/pp(end) (4*s_Boltzman*aa(end)*Kelvin0(end)^3 + alphac(end))/pp(1) ];
    
    for i=1:n+1
        for j=1:n+1
            for k=1:n+1
                for l=1:n+1
                    for o=1:n+1
                        for p=1:n+1
                            p3_par(i,j,k,l,o,p)=(4*s_Boltzman*aa(i)*Kelvin0(o)^3*( 1 - (aa(i)-bb)/aa(i)) ...
                                +alphac(p)+alphae(j)*6.112*exp(Temperatura_rif(k).*17.67/ ...
                                (Temperatura_rif(k)+243.5))*17.67*243.5/ ...
                                (Temperatura_rif(k)+243.5)^2)/pp(l);
                        end
                    end
                end
            end
        end
    end
    p3=[min(min(min(min(min(min(p3_par))))))  max(max(max(max(max(max(p3_par))))))];
    
    p4=[1 100*meanD(prof)^(-0.35)];
    p5=[( (1-rs(end))*minRad_a + alphae(1)*Delta_ea(1)+Delta_alphae(1)*(ea(1) - ew(end) + Delta_ea(1)) + Delta_alphac(1)*DT(1) )/pp(end) ...
        ( (1-rs(1))*maxRad_a + alphae(end)*Delta_ea(end)+Delta_alphae(end)*(ea(end) - ew(1) + Delta_ea(end)) + Delta_alphac(end)*DT(end) )/pp(1)];
    p6=[0 1];
    p7=[0 150];
    p8=[0 0.5];
    
    par(1,prof,:)=p1;
    par(2,prof,:)=p2;
    par(3,prof,:)=p3;
    par(4,prof,:)=p4;
    par(5,prof,:)=p5;
    par(6,prof,:)=p6;
    par(7,prof,:)=p7;
    par(8,prof,:)=p8;
end

par(1,:,2)=min(2,par(1,:,2));
% par(5,:,2)=min(2,par(1,:,2));
par(5,:,1)=max(0,par(1,:,1));

% Regression curves from Toffolon et al., (2014)
par_8(1,:)=0.487944 - 0.096334*log(meanD);
par_4(1,:)=-0.042417 +0.01745*log(meanD);
par_8(2,:)=0.207095*meanD.^(-0.67172);
par_4(2,:)=0.222683*meanD.^(-0.635202);
par_8(3,:)=0.262359*meanD.^(-0.658578);
par_4(3,:)=0.1753*meanD.^(-0.540387);
par_8(4,:)=31.331137*meanD.^(-0.329713);
par_4(4,:)=35.382662*meanD.^(-0.360257);
par_8(5,:)=0.843294*meanD.^(-0.731699);
par_8(6,:)=0.627705 -0.030032* log(meanD);

figure
for i=1:6
    subplot(3,2,i)
    if i==2 | i==3 | i==4 | i==5
        loglog(meanD,par(i,:,1)); hold on
        loglog(meanD,par(i,:,2));
    elseif i==1 | i==6
        semilogx(meanD,par(i,:,1)); hold on
        semilogx(meanD,par(i,:,2));
    end
    grid on
    if i==1
        semilogx(meanD,par_8(1,:),'-k', 'linewidth',2); hold on
        semilogx(meanD,par_4(1,:),'-m', 'linewidth',2); hold on
    elseif i==2
        loglog(meanD,par_8(2,:),'-k', 'linewidth',2); hold on
        loglog(meanD,par_4(2,:),'-m', 'linewidth',2); hold on
    elseif i==3
        loglog(meanD,par_8(3,:),'-k', 'linewidth',2); hold on
        loglog(meanD,par_4(3,:),'-m', 'linewidth',2); hold on
    elseif i==4
        loglog(meanD,par_8(4,:),'-k', 'linewidth',2); hold on
        loglog(meanD,par_4(4,:),'-m', 'linewidth',2); hold on
    elseif i==5
        semilogx(meanD,par_8(5,:),'-k', 'linewidth',2); hold on
    elseif i==6
        semilogx(meanD,par_8(6,:),'-k', 'linewidth',2); hold on
    end
end

set(gcf,'paperunits','centimeters','papersize',[13 20],'paperposition',[0 0 13 20]);
print(gcf,'-r300','-dpdf', 'parameter_range.pdf');
print(gcf,'-r500','-dpng', 'parameter_range.png');


% creo range parametri crescente attorno a valore atteso. Limitato da range
% massimo
Ilake=find(meanD==depth); Ilake=Ilake(1);

par_range=squeeze(par(:,Ilake,:));
P=[par_range(:,1)' ; par_range(:,2)'];
save(['parameters_depth=' num2str(depth) 'm.txt'],'P','-ascii');
