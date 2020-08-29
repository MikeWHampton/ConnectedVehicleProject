function energy = EConsumpFunc8(x0)
alpha01 = x0(1);
beta01 = x0(2);
beta02 = x0(3);
beta03 = x0(4);
beta04 = x0(5);
beta05 = x0(6);
beta06 = x0(7);
gamma01 = x0(8);
gamma02 = x0(7);
gamma03 = x0(10);
gamma04 = x0(11);
gamma05 = x0(12);
gamma06 = x0(13);
% dataset to use
% 4 vehicles, 3 CHV + 1 CAV
% load('20170727_202308_4_vehicles');
% kcav=4;

%8 vehicles, 6 CHV + 1 CAV + 1 CHV
load('20170730_114049_8_vehicles');
kcav=7;

alpha=0.4;
beta=[0.2,0.3,0.3,0,0,0];
% fixed parameters for CAV
sigma=0.6; hst=5; hgo=55; vmax=30; amin=7; amax=3;
betasum=sum(beta);

% parameters for safe and non-conservative bounds
kappa_min=0.35; kappa_max=1.5; hst_min=3; hst_max=10;

% number of measured vehicles
veh_num=length(time);
% vehicle number of CAV: kcav
%alpha01=0.3; beta01=0.4; gamma01=0.3; beta02=0.3; gamma02 = 0.2; % Designed Parameters

% range policy and saturations for CAV
V=@(h)vmax*(hgo<=h) + vmax*(3*hgo-hst-2*h)*((h-hst)^2)/((hgo-hst)^3).*(hst<h & h<hgo);
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;

% control input for CAV
u=@(h0,v0,v1,v2,v3,v4,v5,v6,a1,a2,a3,a4,a5,a6) alpha01*(V(h0)-v0)+beta01*(W(v1)-v0)+...
    beta02*(W(v2)-v0)+beta03*(W(v3)-v0)+beta04*(W(v4)-v0)+beta05*(W(v5)-v0)+beta06*(W(v6)-v0)+...
    gamma01*a1+gamma02*a2+gamma03*a3+gamma04*a4+gamma05*a5+gamma06*a6;
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

v4 = @(t)vLead(t==tsim, 4);
v4delay=@(t)vLeadDelay(t==tsim, 4);

v5 = @(t)vLead(t==tsim, 5);
v5delay = @(t)vLeadDelay(t==tsim, 5);

v6 = @(t)vLead(t==tsim, 6);
v6delay = @(t)vLeadDelay(t==tsim, 6);

a1=@(t)aLead(t==tsim, 1);
a1delay=@(t)aLeadDelay(t==tsim, 1);

a2=@(t)aLead(t==tsim, 2);
a2delay=@(t)aLeadDelay(t==tsim, 2);

a3=@(t)aLead(t==tsim, 3);
a3delay=@(t)aLeadDelay(t==tsim, 3);

a4=@(t)aLead(t==tsim, 4);
a4delay=@(t)aLeadDelay(t==tsim, 4);
    
a5=@(t)aLead(t==tsim, 5);
a5delay=@(t)aLeadDelay(t==tsim, 5);
    
a6=@(t)aLead(t==tsim, 6);
a6delay=@(t)aLeadDelay(t==tsim, 6);

%% Simulation
% right-hand side of equations
model=@(t,x,xdelay)[v1(t)-x(2);
                    sat(u(xdelay(1),xdelay(2),v1delay(t),v2delay(t),v3delay(t),...
                    v4delay(t),v5delay(t),v6delay(t),a1delay(t),a2delay(t),...
                    a3delay(t),a4delay(t),a5delay(t),a6delay(t)))];

% perform simulation
x=ddeab4(@(t,x,xdelay)model(t,x,xdelay),sigma,xinit,tsim);

% extract headway and velocity
hdwysim=x(1,:).';
velsim=x(2,:).';

% check safe and non-conservative bounds
[safe_bound,non_conservative_bound] = check_bounds(hdwysim,velsim,kappa_min,kappa_max,hst_min,hst_max);
if strcmp(safe_bound, 'Safe bound is violated!!')
    energy = 10000;
    return
end
if strcmp(non_conservative_bound, 'Non-conservative bound is violated!!')
    energy = 10000;
    return
end

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

energy = enconssim(end);
end
