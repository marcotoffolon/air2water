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

% POST-PROCESSING
% Aim: to analyze and plot the results

close all
clear all
clc

%% modify the following lines according to your needs
folder='Superior/output_2/';
dt='1d';
IDair='stndrck';
IDwat='satt';
index='RMS';
runmode='PSO';
toll=2;  % minimum efficiency index used to make the dotty plots (maximum if index = RMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..
set(0,'defaulttextfontname','Times New Roman','defaultaxesfontname','Times New Roman',...
    'defaulttextfontsize',10,'defaultaxesfontsize',10);
% Figures are produced using a colorblind barrier-free color pallet (jfly.iam.u-tokyo.ac.jp/color/)
orange=[230 159 0]/255;
blue=[0 114 178]/255;
light_blue=[86 180 233]/255;

%% 1. Dotty plots
file=['0_' runmode '_' index '_' IDair '_' IDwat '_c_' dt '.out'];
fid=fopen([folder file]);
prec='float64';
parset=fread(fid,prec);
parset=reshape(parset,9,length(parset)/9)';
fclose(fid);

eff=parset(:,end);          % efficiency indexes associated to each parameter sets
parset=parset(:,1:end-1);   % parameter sets
if strcmp(index,'RMS')==1
    eff=-eff;   % when RMS is used as efficiency index, air2water works with -RMS
    I=find(eff<=(toll));
    parset=parset(I,:);
    eff=eff(I);
    [best_eff I_best]=min(eff);
    plot_limits=[best_eff*0.9 toll];
else
    I=find(eff>=(toll));
    parset=parset(I,:);
    eff=eff(I);
    [best_eff I_best]=max(eff);
    plot_limits=[toll best_eff*1.1];
end


figure; hold on
for i=1:8
    subplot(2,4,i)
    plot(parset(:,i),eff,'.k'); hold on
    plot(parset(I_best,i),best_eff,'.','color',orange,'markersize',15)
    ylim(plot_limits);
    xlabel(['par' num2str(i)]);
    if strcmp(index,'RMS')==1
        ylabel([index '[°C]'])
    else
        ylabel(index);
    end
end

set(gcf,'paperunits','centimeters','papersize',[18 10],'paperposition',[0 0 18 10]);
print(gcf,'-r300','-dpdf', [folder 'dottyplots_' runmode '_' index '_' IDair '_' IDwat '.pdf']);
print(gcf,'-r300','-dpng', [folder 'dottyplots_' runmode '_' index '_' IDair '_' IDwat '.png']);

%% 2. Series
file_cal=['2_' runmode '_' index '_' IDair '_' IDwat '_cc_' dt '.out'];
file_val=['3_' runmode '_' index '_' IDair '_' IDwat '_cv_' dt '.out'];

T_cal=load([folder file_cal]);
T_cal=T_cal(366:end,:); % the first is a warm up year
T_cal(T_cal==-999)=NaN;
date_cal=datenum([T_cal(:,1:3)]);
RMSE_cal=sqrt(nanmean((T_cal(:,5)-T_cal(:,6)).^2));

figure
plot(date_cal,T_cal(:,4),'.','color',light_blue); hold on
plot(date_cal,T_cal(:,5),'.','color',blue); 
plot(date_cal,T_cal(:,6),'.','color',orange); 
title(['Calibration, RMSE=' num2str(RMSE_cal) '°C'])
xlabel('Time');
ylabel('Temperature [°C]');
legend('Air temperature','Observed water temperature','Simulated water temperature','location','SouthEast')
datetick('x','mmm-yy')
set(gcf,'paperunits','centimeters','papersize',[18 10],'paperposition',[0 0 18 10]);
print(gcf,'-r300','-dpdf', [folder 'calibration_' runmode '_' index '_' IDair '_' IDwat '.pdf']);
print(gcf,'-r300','-dpng', [folder 'calibration_' runmode '_' index '_' IDair '_' IDwat '.png']);

if exist([folder file_val],'file')
    T_val=load([folder file_val]);
    T_val=T_val(366:end,:); % the first is a warm up year
    T_val(T_val==-999)=NaN;
    date_val=datenum([T_val(:,1:3)]);
    RMSE_val=sqrt(nanmean((T_val(:,5)-T_val(:,6)).^2));
    
    figure
    plot(date_val,T_val(:,4),'.','color',light_blue); hold on
    plot(date_val,T_val(:,5),'.','color',blue);
    plot(date_val,T_val(:,6),'.','color',orange);
    title(['Validation, RMSE=' num2str(RMSE_val) '°C'])
    xlabel('Time');
    ylabel('Temperature [°C]');
    legend('Air temperature','Observed water temperature','Simulated water temperature','location','SouthEast')
    datetick('x','mmm-yy')
    set(gcf,'paperunits','centimeters','papersize',[18 10],'paperposition',[0 0 18 10]);
    print(gcf,'-r300','-dpdf', [folder 'validation_' runmode '_' index '_' IDair '_' IDwat '.pdf']);
    print(gcf,'-r300','-dpng', [folder 'validation_' runmode '_' index '_' IDair '_' IDwat '.png']);
end

cd('post_processing')