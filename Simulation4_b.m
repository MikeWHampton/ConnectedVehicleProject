% load('20170727_202308_4_vehicles');
% kcav=4;
% alpha=0.4;
% beta=[0.2,0.3,0.3];
% 8 vehicles, 6 CHV + 1 CAV + 1 CHV
load('20170730_114049_8_vehicles');
kcav=7;
% alpha=0.4;
% beta=[0.2,0.3,0.3,0,0,0];

% fixed parameters for CAV
sigma=0.6; hst=5; hgo=55; vmax=30; amin=7; amax=3;
%betasum=sum(beta);

% parameters for safe and non-conservative bounds
kappa_min=0.35; kappa_max=1.5; hst_min=3; hst_max=10;

% number of measured vehicles
veh_num=length(time);
% vehicle number of CAV: kcav
alpha01=  0.016012698808363; 
beta01=   0.0605997101698718; 
beta02=   0.0790348238573579; 
beta03 =  0.0720348241593848; 
gamma01=  -2.03124251184872e-05; 
gamma02 = 0.00286837257041867; 
gamma03 = 0.204701500130322; 
gamma04 = 0.104549966187748; 
gamma05 = 0.0514667972784925; 
gamma06 = 0.241430413987001;% Designed Parameters

% range policy and saturations for CAV
V=@(h)vmax*(hgo<=h) + vmax.*(3*hgo-hst-2.*h)*((h-hst).^2)/((hgo-hst).^3).*(hst<h & h<hgo);
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;

% control input for CAV
u=@(h0,v0,v1,v2,v3,a1,a2,a3,a4,a5,a6)alpha01*(V(h0)-v0)+beta01*(W(v1)-v0)+...
    beta02*(W(v2)-v0)+beta03*(W(v3)-v0)+gamma01*a1+gamma02*a2+gamma03*a3+...
    gamma04*a4+gamma05*a5+gamma06*a6;
% initial conditions for CAV
h0=hdwy{kcav}(1);
v0=vel{kcav}(1);
xinit=@(t)[h0;v0];

% simulation time
t0=min(vertcat(time{:}));
tend=max(vertcat(time{:}));
deltat=min(diff(time{1}));
tsim=(t0:deltat:tend).';        

% number of preceding CHVs that the CAV may respond to
% number of preceding CHVs that the CAV may respond to
chv_num=length(beta);
% preceding CHVs' velocity from data
for kL=1:chv_num
    vLead(:,kL)=interp1(time{kcav-kL},vel{kcav-kL},tsim,'linear','extrap');
    vLeadDelay(:,kL)=interp1(tsim,vLead(:,kL),tsim-sigma,'linear','extrap');
    aLead(:, kL)=interp1(time{kcav-kL},acc{kcav-kL},tsim,'linear','extrap');
    aLeadDelay(:, kL) =interp1(tsim,aLead(:,kL),tsim-sigma,'linear','extrap');
end

v1 = @(t)vLead(t==tsim, 1);
v1delay=@(t)vLeadDelay(t==tsim, 1);

v2 = @(t)vLead(t==tsim, 2);
v2delay = @(t)vLeadDelay(t==tsim, 2);

v3 = @(t)vLead(t==tsim, 3);
v3delay = @(t)vLeadDelay(t==tsim, 3);

a1=@(t)aLead(t==tsim, 1);
a1delay=@(t)aLeadDelay(t==tsim, 1);

a2=@(t)aLead(t==tsim, 2);
a2delay=@(t)aLeadDelay(t==tsim, 2);

a3=@(t)aLead(t==tsim, 3);
a3delay=@(t)aLeadDelay(t==tsim, 3);

if kcav == 7
    a4=@(t)aLead(t==tsim, 4);
    a4delay=@(t)aLeadDelay(t==tsim, 4);
    
    a5=@(t)aLead(t==tsim, 5);
    a5delay=@(t)aLeadDelay(t==tsim, 5);
    
    a6=@(t)aLead(t==tsim, 6);
    a6delay=@(t)aLeadDelay(t==tsim, 6);
end

% legend to put on figure
vehtype=repmat({'measured CHV'},1,veh_num);
vehtype{kcav}='measured CAV';
vehlegend=horzcat(vehtype,{'simulated CAV'});

%% Simulation of MY Controller
% right-hand side of equations
model=@(t,x,xdelay)[v1(t)-x(2);
                    sat(u(xdelay(1),xdelay(2),v1delay(t),v2delay(t),...
                    v3delay(t),a1delay(t),a2delay(t),a3delay(t),a4delay(t),a5delay(t),a6delay(t)))];

% perform simulation
x=ddeab4(@(t,x,xdelay)model(t,x,xdelay),sigma,xinit,tsim);

% extract headway and velocity
hdwysim1=x(1,:).';
velsim1=x(2,:).';
%% Initialize THEIR Controller
kcav=2; alpha01=0.4; beta01=0.2; beta02=0.3; beta03=0.3; K0 = 0.6;
V=@(h)vmax*(hgo<=h) + K0*(h-hst).*(hst<h & h<hgo)
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;

% control input for CAV
u=@(h,v,v1,v2,v3)alpha01*(V(h)-v)+beta01*(W(v1)-v)+beta02*(W(v2)-v)+beta03*(W(v3)-v);

model=@(t,x,xdelay)[v1(t)-x(2);
                    sat(u(xdelay(1),xdelay(2),v1delay(t),v2delay(t),v3delay(t)))];
% perform simulation
x2=ddeab4(@(t,x,xdelay)model(t,x,xdelay),sigma,xinit,tsim);
hdwysim2=x2(1,:).';
velsim2=x2(2,:).';

% check safe and non-conservative bounds
check_bounds(hdwysim1,velsim1,kappa_min,kappa_max,hst_min,hst_max);

% calculate accelearation
accsim1=sgolayfilt(diff(velsim1)/deltat,3,21);
accsim1=[accsim1(1);accsim1];
accsim2=sgolayfilt(diff(velsim2)/deltat,3,21);
accsim2=[accsim2(1);accsim2];
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
enconssim1=cumsum(max(accsim1+a+c*velsim1.^2,0).*velsim1*deltat);
enconssim2=cumsum(max(accsim2+a+c*velsim2.^2,0).*velsim2*deltat);
%% Plot of solution
% plot velocity vs time
figure(1); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot speed of simulated CAV
LL(end)=plot(tsim,velsim1,'Linewidth',1.5,'Color','r');
LL(end)=plot(tsim,velsim2,'Linewidth',1.5,'Color','b');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Speed [m/s]');
title('Speed vs Time');
legend('designed controller','baseline controller','Location','northwest');

% plot acceleration vs time
figure(2); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot acceleration of simulated vehicle
LL(end)=plot(tsim,accsim1,'Linewidth',1.5,'Color','r');
LL(end)=plot(tsim,accsim2,'Linewidth',1.5,'Color','b');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Acceleration [m/s^2]');
legend('designed controller','baseline controller','Location','northwest');

% plot headway vs time
figure(3); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot headway of simulated vehicle
LL(end)=plot(tsim,hdwysim1,'Linewidth',1.5,'Color','r');
LL(end)=plot(tsim,hdwysim2,'Linewidth',1.5,'Color','b');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Headway [m]');
legend('designed controller','baseline controller','Location','northwest');

% plot energy consumption vs time
figure(4); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot energy consumption of simulated vehicle
LL(end)=plot(tsim,enconssim1,'Linewidth',1.5,'Color','r');
LL(end)=plot(tsim,enconssim2,'Linewidth',1.5,'Color','b');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Energy consumption [J/kg]');
legend('designed controller','baseline controller','Location','northwest');

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
velFFT1=cell(size(vel));
velFFT2=cell(size(vel));
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
velinterp1=interp1(tsim,velsim1,tFFT,'linear','extrap');
velsimFFT1=abs(fft(velinterp1-mean(velinterp1))/n);
velsimFFT1=velsimFFT1(1:floor(n/2+1));
velsimFFT1=2*velsimFFT1;
velsimFFT1=sgolayfilt(velsimFFT1,3,31);

velinterp2=interp1(tsim,velsim2,tFFT,'linear','extrap');
velsimFFT2=abs(fft(velinterp2-mean(velinterp2))/n);
velsimFFT2=velsimFFT2(1:floor(n/2+1));
velsimFFT2=2*velsimFFT2;
velsimFFT2=sgolayfilt(velsimFFT2,3,31);
% plot frequency content
figure(5); clf; box on;
semilogy(freq,velsimFFT1,'Linewidth',1.5,'Color','r');
hold on
semilogy(freq,velsimFFT2,'Linewidth',1.5,'Color','b');
xlim([0,1]);
% ylim([0.001,1]);
xlabel('Frequency [Hz]');ylabel('Velocity [m/s]');
legend('designed controller','baseline controller','Location','northwest');

% plot frequency response
figure(6); clf; box on;
semilogy(freq,velsimFFT1./velFFT{1},'Linewidth',1.5,'Color','r');
hold on
semilogy(freq,velsimFFT2./velFFT{1},'Linewidth',1.5,'Color','b');
xlim([0,1]);
% ylim([0.1,10]);
xlabel('Frequency [Hz]');ylabel('Amplification');
legend('designed controller','baseline controller','Location','northwest');

% string stability index
Cs1=sum(max(velsimFFT1./velFFT{1}-1,0)*fs/n);
disp('String stability index for designed controller: ')
disp(Cs1)
Cs2=sum(max(velsimFFT2./velFFT{1}-1,0)*fs/n);
disp('String stability index for baseline controller: ')
disp(Cs2)

%% Evaluation of time to collision
% calculate time to collision
ttc=cell(size(vel));
for kk=2:veh_num
    ttc{kk}=interp1(tCom{kk},hdwy{kk},time{kk})./...
            (vel{kk}-interp1(time{kk-1},vel{kk-1},time{kk}));
    ttc{kk}(ttc{kk}<0)=nan;
end
ttcsim1=hdwysim1./(velsim1-interp1(time{kcav-1},vel{kcav-1},tsim));
ttcsim2=hdwysim2./(velsim2-interp1(time{kcav-1},vel{kcav-1},tsim));
ttcsim1(ttcsim1<0)=nan;
ttcsim2(ttcsim2<0)=nan;
% plot time to collision vs time
figure(7); clf; hold on; box on;
plot(tsim,ttcsim1,'Linewidth',1.5,'Color','r');
hold on
plot(tsim,ttcsim2,'Linewidth',1.5,'Color','b');
xlim([t0,tend]);
ylim([0,20]);
xlabel('Time [s]');ylabel('Time to collision [s]');
legend('designed controller','baseline controller','Location','northwest');

% time to collision index
ttccr=2;
Ct=sum(max(ttccr-ttcsim1,0)*deltat);

%% Evaluation of headway error
% calculate headway error (difference from desired headway)
Vinv=@(v)hst+v/vmax*(hgo-hst);  % valid for 0<v<vmax only
hdwyerr=cell(size(vel));
for kk=2:veh_num
    hdwyerr{kk}=hdwy{kk}-interp1(time{kk},Vinv(vel{kk}),tCom{kk});
end
hdwyerrsim1=hdwysim1-Vinv(velsim1);
hdwyerrsim2=hdwysim2-Vinv(velsim2);

% plot headway error vs time
figure(8); clf; hold on; box on;
plot(tsim,hdwyerrsim1,'Linewidth',1.5,'Color','r');
plot(tsim,hdwyerrsim2,'Linewidth',1.5,'Color','b');
xlim([t0,tend]);
xlabel('Time [s]');ylabel('Headway error [m]');
legend('designed controller','baseline controller','Location','northwest');

%% Plot bounds and data on (h,v) plane for the simulated controller
figure(9); clf; hold on; box on;
plot(hdwysim1,velsim1,'Linewidth',1.5,'Color','r')
hold on
plot(hdwysim2,velsim2,'Linewidth',1.5,'Color','b')
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
legend('designed controller','baseline controller','Location','northwest');
