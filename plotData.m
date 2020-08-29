% Main data:
% veh_num - number of vehicles
% time - time
% vel - velocity
% acc - acceleration
% arc - arclength (starts from zero for the leading car)
% encons - energy consumption
% tCom - common time instants when data is available for both vehicles in a pair
% hdwy - headway from predecessor vehicle at these common time instants
%        (great circle distance - 5 meters vehicle length)
% ...Loss - data filled with NaNs at packet losses (diff(time)>0.2 or diff(tCom)>0.2 sec)
%           this way, the plotted line is interrupted at packet losses
% ...Rec - last data received before packet losses
%          and first data received after packet losses
% time0 - measurement start time
% tstart, tend - start and end time for the relevant part of data being plotted

clear; clc; %close all;

% Select data set to be loaded
% load('20170611_210016_2_vehicles.mat');
% load('20170727_202308_4_vehicles.mat');
load('20170730_114049_8_vehicles.mat');

% Number of vehicles
veh_num=size(colours,1);

% Plot arclength vs time
figure(1);clf;hold on;box on;grid on;
for kk=1:veh_num
%    plot(time{kk},arc{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeLoss{kk},arcLoss{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeRec{kk},arcRec{kk},'Linewidth',1.5,'Color',colours(kk,:),'LineStyle','none','Marker','x');
end
xlabel('Time [s]');ylabel('Arclength [m]');
xlim([tstart tend]);

% Plot velocity vs time
figure(2);clf;hold on;box on;grid on;
for kk=1:veh_num
%    plot(time{kk},vel{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeLoss{kk},velLoss{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeRec{kk},velRec{kk},'Linewidth',1.5,'Color',colours(kk,:),'LineStyle','none','Marker','x');
end
xlabel('Time [s]');ylabel('Speed [m/s]');
xlim([tstart tend]);

% Plot acceleration vs time
figure(3);clf;hold on;box on;grid on;
for kk=1:veh_num
%    plot(time{kk},acc{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeLoss{kk},accLoss{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeRec{kk},accRec{kk},'Linewidth',1.5,'Color',colours(kk,:),'LineStyle','none','Marker','x');
end
xlabel('Time [s]');ylabel('Acceleration [m/s^2]');
xlim([tstart tend]);

% Plot headway vs time
figure(4);clf;hold on;box on;grid on;
for kk=1:veh_num
%    plot(tCom{kk},hdwy{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(tComLoss{kk},hdwyLoss{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(tComRec{kk},hdwyRec{kk},'Linewidth',1.5,'Color',colours(kk,:),'LineStyle','none','Marker','x');
end
xlabel('Time [s]');ylabel('Headway [m]');
xlim([tstart tend]);

% Plot energy consumption vs time
figure(5);clf;hold on;box on;grid on;
for kk=1:veh_num
%    plot(time{kk},encons{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeLoss{kk},enconsLoss{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeRec{kk},enconsRec{kk},'Linewidth',1.5,'Color',colours(kk,:),'LineStyle','none','Marker','x');
end
xlabel('Time [s]');ylabel('Energy consumption [J/kg]');
xlim([tstart tend]);

% Plot bounds and data on (h,v) plane for the measured CAV
figure(6);clf;hold on;box on;grid on;

% parameters for safe and non-conservative bounds
kappa_min=0.35;
kappa_max=1.5;
hst_min=3;
hst_max=10;

if length(time) == 8
    hdwy_CAV = interp1(tCom{end-1},hdwy{end-1},time{end-1});
    plot(hdwy_CAV,vel{end-1},'Linewidth',2)
else
    hdwy_CAV = interp1(tCom{end},hdwy{end},time{end});
    plot(hdwy_CAV,vel{end},'Linewidth',2)
end
hold on
h = -10:0.01:100;
vsafe = (h-hst_min)*kappa_max;
vncon = (h-hst_max)*kappa_min;
plot(h,vsafe,'r','Linewidth',2)
hold on
plot(h,vncon,'r','Linewidth',2)
xlim([0 70])
ylim([0 30])
xlabel('Headway h','Fontsize',16)
ylabel('Velocity v','Fontsize',16)
title(['Measured CAV in ',num2str(veh_num),'-vehicle dataset'],'Fontsize',16)