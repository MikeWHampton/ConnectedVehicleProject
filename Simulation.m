% Perform simulations for a 2 vehicle chain using data as input
clear; clc;
close all;

%% Parameters
load('20170611_210016_2_vehicles');
kcav=2; alpha01=0.1; beta01=0.25; gamma01=0.3; % Designed Parameters
% fixed parameters for CAV
sigma=0.6; hst=5; hgo=55; vmax=30; amin=7; amax=3;

% parameters for safe and non-conservative bounds
kappa_min=0.35; kappa_max=1.5; hst_min=3; hst_max=10;

% number of measured vehicles
veh_num=length(time);
% vehicle number of CAV: kcav

% range policy and saturations for CAV
V=@(h)vmax*(hgo<=h) + vmax*(3*hgo-hst-2*h)*((h-hst)^2)/((hgo-hst)^3).*(hst<h & h<hgo);
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;

% control input for CAV
%OVM with acceleration feedback
u=@(h,v,vL,aL, hdot)alpha01*(V(h)-v)+beta01*(W(vL)-v)+gamma01*aL;
%IDM
%Overview: http://www-personal.umich.edu/~orosz/articles/ITS_2017_Jin.pdf
%More in depth with parameters: https://arxiv.org/pdf/cond-mat/0002177.pdf
%a = 2; b = 1.5; Tgap = 1;
%u=@(h,v,vL,aL, hdot)a*(1-(v/vmax)^4-((hst+Tgap*v-(hdot)*v/sqrt(4*a*b))/h)^2);

% initial conditions for CAV
h0=hdwy{kcav}(1);
v0=vel{kcav}(1);
xinit=@(t)[h0;v0];

% simulation time
t0=min(vertcat(time{:}));
tend=max(vertcat(time{:}));
deltat=min(diff(time{1}));
tsim=(t0:deltat:tend).';

% leader's velocity from data
vLead=interp1(time{kcav-1},vel{kcav-1},tsim,'linear','extrap');
vLeaddelay=interp1(tsim,vLead,tsim-sigma,'linear','extrap');
vL=@(t)vLead(t==tsim);
vLdelay=@(t)vLeaddelay(t==tsim);
aLead=interp1(time{kcav-1},acc{kcav-1},tsim,'linear','extrap');
aLeaddelay=interp1(tsim,aLead,tsim-sigma,'linear','extrap');
aL=@(t)aLead(t==tsim);
aLdelay=@(t)aLeaddelay(t==tsim);

% title to put on figure
problem='CCC using measurement data';
% list of parameters to put on figure
parlist=['\alpha=',num2str(alpha01,'%3.2f'),' [1/s]   ',...
         '\beta=',num2str(beta01,'%3.2f'),' [1/s]   ',...
         '\gamma_0_1=',num2str(gamma01,'%3.2f'),' [1/s^2]   ',...
         '\sigma=',num2str(sigma,'%3.2f'),' [s]'];
% parlist=['\alpha=',num2str(a,'%3.2f'),'    ',...
%          '\beta=',num2str(b,'%3.2f'),'    ',...
%          'T_g_a_p=',num2str(Tgap,'%3.2f'),'    ',...
%          '\sigma=',num2str(sigma,'%3.2f')];

% legend to put on figure
vehtype=repmat({'measured CHV'},1,veh_num);
vehtype{kcav}='measured CAV';
vehlegend=horzcat(vehtype,{'simulated CAV'});

%% Simulation
% right-hand side of equations
model=@(t,x,xdelay)[vL(t)-x(2);
                    sat( u(xdelay(1), xdelay(2), vLdelay(t), aLdelay(t), vLdelay(t)-xdelay(2))  )];

% perform simulation
x=ddeab4(@(t,x,xdelay)model(t,x,xdelay),sigma,xinit,tsim);

% extract headway and velocity
hdwysim=x(1,:).';
velsim=x(2,:).';

% check safe and non-conservative bounds
check_bounds(hdwysim,velsim,kappa_min,kappa_max,hst_min,hst_max);

% calculate accelearation
accsim=sgolayfilt(diff(velsim)/deltat,3,21);
accsim=[accsim(1);accsim];
    
% calculate energy consumption
gamma=0.01;         % [-] tyre rolling resistance coefficient
g=9.81;             % [m/s^2] gravitatioinal constant
a=gamma*g;          % [m/s^2]
Cd=0.34;            % [-] air drag coefficient
A=2.32;             % [m^2] frontal area
rho=1.23;           % [kg/m^3] air density at 25 degree
k=0.5*Cd*rho*A;     % [kg/m]
m=1770;             % [kg] mass of the vehicle
c=k/m;              % [1/m]
enconssim=cumsum(max(accsim+a+c*velsim.^2,0).*velsim*deltat);

%% Plot of solution
% plot velocity vs time
figure(1); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot measured velocity
for kk=1:veh_num
%    plot(time{kk},vel{kk},'Linewidth',1.5,'Color',colours(kk,:));
   LL(kk)=plot(timeLoss{kk},velLoss{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeRec{kk},velRec{kk},'Linewidth',1.5,'Color',colours(kk,:),'LineStyle','none','Marker','x');
end
% plot speed of simulated CAV
LL(end)=plot(tsim,velsim,'Linewidth',2.5,'Color','r');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Speed [m/s]');
title(['Simulations of ',problem,' - ',parlist]);
legend(LL,vehlegend,'Location','northwest');

% plot acceleration vs time
figure(2); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot measured acceleration
for kk=1:veh_num
%    plot(time{kk},acc{kk},'Linewidth',1.5,'Color',colours(kk,:));
   LL(kk)=plot(timeLoss{kk},accLoss{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeRec{kk},accRec{kk},'Linewidth',1.5,'Color',colours(kk,:),'LineStyle','none','Marker','x');
end
% plot acceleration of simulated vehicle
LL(end)=plot(tsim,accsim,'Linewidth',2.5,'Color','r');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Acceleration [m/s^2]');
legend(LL,vehlegend,'Location','northwest');

% plot headway vs time
figure(3); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot measured headway
for kk=1:veh_num
%    plot(tCom{kk},hdwy{kk},'Linewidth',1.5,'Color',colours(kk,:));
   LL(kk)=plot(tComLoss{kk},hdwyLoss{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(tComRec{kk},hdwyRec{kk},'Linewidth',1.5,'Color',colours(kk,:),'LineStyle','none','Marker','x');
end
% plot headway of simulated vehicle
LL(end)=plot(tsim,hdwysim,'Linewidth',2.5,'Color','r');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Headway [m]');
legend(LL,vehlegend,'Location','northwest');

% plot energy consumption vs time
figure(4); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot measured energy consumption
for kk=1:veh_num
%    plot(time{kk},encons{kk},'Linewidth',1.5,'Color',colours(kk,:));
   LL(kk)=plot(timeLoss{kk},enconsLoss{kk},'Linewidth',1.5,'Color',colours(kk,:));
   plot(timeRec{kk},enconsRec{kk},'Linewidth',1.5,'Color',colours(kk,:),'LineStyle','none','Marker','x');
end
% plot energy consumption of simulated vehicle
LL(end)=plot(tsim,enconssim,'Linewidth',2.5,'Color','r');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Energy consumption [J/kg]');
legend(LL,vehlegend,'Location','northwest');

%% Evaluation of frequency response
% for FFT, consider only those times when vehicles were moving
% when the CAV's speed (and possibly its predecessor's speed) exceeded vmin
vmin=1;
% get time interval (t1<=t<=t2) corresponding to moving
t1=time{kcav}(find(vel{kcav}>=vmin & ...
    interp1(time{kcav-1},vel{kcav-1},time{kcav})>=vmin,1,'first'));
t2=time{kcav}(find(vel{kcav}>=vmin,1,'last'));

% pick time vector to which data is interpolated before doing FFT
tFFT=tsim(t1<=tsim & tsim<=t2);
% get number of data points
n=length(tFFT);
% get sampling frequency and discrete frequency values
fs=1/deltat;
freq=(0:floor(n/2))/n*fs;

% frequency content of measured velocity
velFFT=cell(size(vel));
for kk=1:veh_num
    % interpolate velocities to the same time vector
    velinterp=interp1(time{kk},vel{kk},tFFT,'linear','extrap');
    % calculate FFT for velocity fluctuations
    velFFT{kk}=abs(fft(velinterp-mean(velinterp))/n);
    % take twice the first half of the FFT
    velFFT{kk}=velFFT{kk}(1:floor(n/2+1));
    velFFT{kk}(2:end)=2*velFFT{kk}(2:end);
    % filter the FFT for smoother results
    velFFT{kk}=sgolayfilt(velFFT{kk},3,31);
end

% frequency content of simulated velocity
velinterp=interp1(tsim,velsim,tFFT,'linear','extrap');
velsimFFT=abs(fft(velinterp-mean(velinterp))/n);
velsimFFT=velsimFFT(1:floor(n/2+1));
velsimFFT(2:end)=2*velsimFFT(2:end);
velsimFFT=sgolayfilt(velsimFFT,3,31);

% plot frequency content
figure(5); clf; box on;
for kk=1:veh_num
    semilogy(freq,velFFT{kk},'Linewidth',1.5,'Color',colours(kk,:));
    hold on;
end
semilogy(freq,velsimFFT,'Linewidth',2.5,'Color','r');
xlim([0,1]);
% ylim([0.001,1]);
xlabel('Frequency [Hz]');ylabel('Velocity [m/s]');
legend(vehlegend,'Location','northeast');

% plot frequency response
figure(6); clf; box on;
for kk=2:veh_num
    semilogy(freq,velFFT{kk}./velFFT{1},'Linewidth',1.5,'Color',colours(kk,:));
    hold on;
end
semilogy(freq,velsimFFT./velFFT{1},'Linewidth',2.5,'Color','r');
xlim([0,1]);
% ylim([0.1,10]);
xlabel('Frequency [Hz]');ylabel('Amplification');
legend(vehlegend(2:end),'Location','northeast');

% string stability index
Cs=sum(max(velsimFFT./velFFT{1}-1,0)*fs/n);
disp('String stability index: ')
disp(Cs)

%% Evaluation of time to collision
% calculate time to collision
ttc=cell(size(vel));
for kk=2:veh_num
    ttc{kk}=interp1(tCom{kk},hdwy{kk},time{kk})./...
            (vel{kk}-interp1(time{kk-1},vel{kk-1},time{kk}));
    ttc{kk}(ttc{kk}<0)=nan;
end
ttcsim=hdwysim./(velsim-interp1(time{kcav-1},vel{kcav-1},tsim));
ttcsim(ttcsim<0)=nan;

% plot time to collision vs time
figure(7); clf; hold on; box on;
for kk=2:veh_num
   plot(time{kk},ttc{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
plot(tsim,ttcsim,'Linewidth',2.5,'Color','r');
xlim([t0,tend]);
ylim([0,20]);
xlabel('Time [s]');ylabel('Time to collision [s]');
legend(vehlegend(2:end),'Location','northeast');

% time to collision index
ttccr=2;
Ct=sum(max(ttccr-ttcsim,0)*deltat);

%% Evaluation of headway error
% calculate headway error (difference from desired headway)
Vinv=@(v)hst+v/vmax*(hgo-hst);  % valid for 0<v<vmax only
hdwyerr=cell(size(vel));
for kk=2:veh_num
    hdwyerr{kk}=hdwy{kk}-interp1(time{kk},Vinv(vel{kk}),tCom{kk});
end
hdwyerrsim=hdwysim-Vinv(velsim);

% plot headway error vs time
figure(8); clf; hold on; box on;
for kk=2:veh_num
   plot(tCom{kk},hdwyerr{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
plot(tsim,hdwyerrsim,'Linewidth',2.5,'Color','r');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Headway error [m]');
legend(vehlegend(2:end),'Location','northeast');

%% Plot bounds and data on (h,v) plane for the simulated controller
figure(9); clf; hold on; box on;
plot(hdwysim,velsim,'Linewidth',2)
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
title(['Simulated controller for ',num2str(veh_num),'-vehicle dataset'],'Fontsize',16)