% INSTRUCTIONS:
% 
% Run every section of the code separately, in the given order
% Display will say which simulation to run and when
%     example: "execute detumble", run project_detumble.slx
% Once the simulation is completed, advance to next section to visualize graphs
% Once graphs have been visualized, advance to next section
% 
% When running project_pointing.slx, first run the sim with pointing control
% toggled off to simulate idle control. Visualize graphs, then run the sim
% again with pointing control toggled on.
% 
% NOTE: project_pointing_RW.slx is to be run instead of project_pointing.slx
% in the sequence of execution is reaction wheel performance is to be 
% evaluated

clear
close
clc

% ORBITAL PARAMETERS

% minimum orbital height [km]
H=515;
% orbital eccentricity [-]
e=0;
% orbital plane inclination [rad]
i=deg2rad(97.5);
% right ascension of ascending node [rad]
OM=deg2rad(0);
% anomaly of pericenter [rad]
om=deg2rad(0); 

%--------------------------------------------------------------------------

mu=astroConstants(13);  % Earth gravitational parameter [km^3/s^2]
rt=astroConstants(23);  % mean Earth radius [km]
J2=1.082e-3;            % coefficient for the second zonal term
r_p=rt+H;               % distance at pericenter [km]
a=r_p/(1-e);            % semi-major axis [km]
n=sqrt(mu/(a^3));       % orbit mean angular velocity [rad/s]
r_s=astroConstants(2);  % mean Sun distance 1AU [km]
th_e=acos(rt/r_s);      % maximum illuminated longitude [rad]
th_a=asin(rt/r_s);      % Sun radiation cone semi-aperture [rad]
T_o=2*pi/n;             % S/C orbital period [s]
T_y=365*24*60*60;       % solar year [s]
T_e=24*60*60;           % Earth day [s]
n_sun=2*pi/T_y;         % Sun mean angular velocity [rad/s]
omega_e=2*pi/T_e;       % Earth rotation angular velocity [rad/s]
eps=deg2rad(23.45);     % ecliptic angle [rad]
m=[0.37 0.32 0.28];     % S/C magnetic dipole [Am^2]

%--------------------------------------------------------------------------

% DYNAMIC PROPERTIES

% dimension along x-axis [m]
dx=0.8;
% dimension along x-axis [m]
dy=0.7;
% dimension along x-axis [m]
dz=0.6;
% mass [kg]
mSC=94;

%--------------------------------------------------------------------------

Ix=mSC/12*(dy^2+dz^2);  % x-axis inertia [kg*m^2]
Iy=mSC/12*(dx^2+dz^2);  % y-axis inertia [kg*m^2]
Iz=mSC/12*(dx^2+dy^2);  % z-axis inertia [kg*m^2]
J=diag([Ix,Iy,Iz]);     % inertial matrix [kg*m^2]
Ky=(Iz-Iy)/Ix;          % inertial yaw coefficients [-]
Kr=(Iz-Ix)/Iy;          % inertial roll coefficients [-]
Kp=(Iy-Ix)/Iz;          % inertial pitch coefficients [-]

%--------------------------------------------------------------------------

% STAR SENSORS

FOV=deg2rad(21.7);              % sensor field of view [rad]
sensor_noise=deg2rad(50/3600);  % maximum sensor noise [rad] 
sensor_sample_time=1;           % sensor sample time [s]

%--------------------------------------------------------------------------

% Rotational Matrices

Rx=@(a) [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
Ry=@(a) [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];
Rz=@(a) [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
   
% x-side sensor installation error
A_sb_1=[999.999999992226e-003    3.46376307675883e-006   -1.88420189887425e-006;
   -3.46376307675883e-006    999.999999994001e-003    3.26316751397826e-012;
    1.88420189887425e-006    3.26316751397826e-012    999.999999998225e-003];
     
% y-side sensor installation error
A_sb_2=[999.999999948501e-003   -10.1488410044293e-006   -22.0178875132149e-012;
    10.1488410044293e-006    999.999999939087e-003    4.33899949057153e-006;
   -22.0178875132149e-012   -4.33899949057153e-006    999.999999990587e-003];
       
% z-side sensor installation error
A_sb_3=[999.999999975613e-003   -46.6812699606578e-012    6.98384920637069e-006;
   -46.6812699606578e-012    999.999999910644e-003    13.3683438779566e-006;
   -6.98384920637069e-006   -13.3683438779566e-006    999.999999886257e-003];

% Star Collection

o=2;    % o=0 12 stars; o=1 42 stars; o=2 162 stars
[S,f]=icosphere(o);
S=S';

% % Discrete Filter
% 
% Ts=sensor_sample_time;  % Sensor sample time [s]
% fs=1/Ts;                % Sensor frequency [Hz]
% ws=2*pi*fs;             % Sensor frequency [rad/s]
% wc=0.01*ws;             % Filter cut frequency [rad/s]
% 
% alpha_f=Ts/(Ts+1/wc);
% Afd=diag((1-alpha_f)*ones(9,1));
% Bfd=Afd;
% Cfd=diag(alpha_f*ones(9,1));
% Dfd=Cfd;

%--------------------------------------------------------------------------

% Magnetic Model

% G matrix
g1_0=-29404.8*1e-9;
g1_1=-1450.9*1e-9;
g2_0=-2499.6*1e-9;
g2_1=2982.0*1e-9;
g2_2=1672*1e-9;
g3_0=1363.2*1e-9;
g3_1=-2381.2*1e-9;
g3_2=1236.2*1e-9;
g3_3=525.7*1e-9;
g4_0=903*1e-9;
g4_1=809.5*1e-9;
g4_2=86.3*1e-9;
g4_3=-309.4*1e-9;
g4_4=48*1e-9;
g5_0=-234.3*1e-9;
g5_1=363.2*1e-9;
g5_2=187.8*1e-9;
g5_3=-140.7*1e-9;
g5_4=-151.2*1e-9;
g5_5=13.5*1e-9;

g_coefficients=[g1_0 g1_1 0     0    0    0;
                g2_0 g2_1 g2_2  0    0    0;
                g3_0 g3_1 g3_2  g3_3 0    0;
                g4_0 g4_1 g4_2  g4_3 g4_4 0;
                g5_0 g5_1 g5_2  g5_3 g5_4 g5_5];

% H matrix
h1_0=0;
h1_1=4652.5*1e-9;
h2_0=0*1e-9;
h2_1=-2991.6*1e-9;
h2_2=-734.6*1e-9;
h3_0=0;
h3_1=-82.1*1e-9;
h3_2=241.9*1e-9;
h3_3=-543.4*1e-9;
h4_0=0;
h4_1=281.9*1e-9;
h4_2=-158.4*1e-9;
h4_3=199.7*1e-9;
h4_4=-349.7*1e-9;
h5_0=0;
h5_1=47.7*1e-9;
h5_2=208.3*1e-9;
h5_3=-121.2*1e-9;
h5_4=32.3*1e-9;
h5_5=98.9*1e-9;

h_coefficients=[h1_0 h1_1 0     0    0    0;
                h2_0 h2_1 h2_2  0    0    0;
                h3_0 h3_1 h3_2  h3_3 0    0;
                h4_0 h4_1 h4_2  h4_3 h4_4 0;
                h5_0 h5_1 h5_2  h5_3 h5_4 h5_5];

% S matrix
S0=1;
S_matrix(1,1)=S0;

for n_index=2:5
    S_matrix(n_index,1)=S_matrix(n_index-1,1)*(2*n_index-1)/n_index;
end

for n_index=1:5
    for m_index=2:6
        if m_index==2
            S_matrix(n_index,m_index)=S_matrix(n_index,m_index-1)*((2)*(n_index-(m_index-1)+1)/ ...
                (n_index+m_index-1))^0.5;
        else
            S_matrix(n_index,m_index)=S_matrix(n_index,m_index-1)*((1)*(n_index-(m_index-1)+1)/ ...
                (n_index+m_index-1))^0.5;
        end
    end
end

g_nm=S_matrix.*g_coefficients;
h_nm=S_matrix.*h_coefficients;

% K matrix
K_matrix=zeros(5,6);
for n_index=2:5
    for m_index=0:n_index
        K_matrix(n_index,m_index+1)=((n_index-1)^2-m_index^2)/((2*n_index-1)*(2*n_index-3));
    end
end

P00=1;
alphaG_0=pi/1.5;    % initial right ascension of Greenwich meridian [rad]

%--------------------------------------------------------------------------

% DETUMBLE CONTROL PARAMETERS

% proportional gain for w [-]
KP=0.035;
% derivative gain for w [-]
KD=0.004;

% SLEW CONTROL PARAMETERS

% proportional gain for alpha_x [-]
kpx=0.6*0.001;
% proportional gain for alpha_y [-]
kpy=kpx;
% proportional gain for alpha_z [-]
kpz=kpy;
% derivative gain for alpha_x [-]
kdx=0.001*450/8;
% derivative gain for alpha_y [-]
kdy=kdx;
% derivative gain for alpha_z [-]
kdz=kdy;

% POINTING CONTROL

% pole placement for state observer
p=-[0.04 0.04 0.04 0.05 0.05+0.025*1i 0.05-0.025*1i];
% nadir pointing angle precision [deg]
nad=2;
% nadir pointing yaw limit [deg]
yaw=20;
% nadir pointing control limit [Nm]
con=1e-4;

%--------------------------------------------------------------------------

% gain control matrices
q=[deg2rad(yaw) deg2rad(1.9*nad) deg2rad(1.9*nad)];
Qc=inv(diag(q))^2;
Rc=(1/(con)^2)*eye(3);

% linearized system matrices
Ac=[0 (1-Ky)*n 0 -Ky*n^2 0 0;
    (Kr-1)*n 0 0 0 -4*Kr*n^2 0;
    0 0 0 0 0 -3*Kp*n^2;
    1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0];

Bc=[J^-1; zeros(3,3)];

Cc=[0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];

Dc=zeros(3,3);

SC=ss(Ac,Bc,Cc,Dc);

Obs=[Bc Ac*Bc (Ac^2)*Bc];            % Observability Matrix
Con=[Cc' Ac'*Cc' (Ac'^2)*Cc'];       % Controllability Matrix

K=lqry(SC,Qc,Rc);       % proportional control gain matrix

Lt=place(Ac',Cc',p);
L=Lt';                  % observer matrix

%--------------------------------------------------------------------------

% THRUSTERS

% thrusters configuration matrix
R=[0 0 dz/2 0 0 -dz/2;
   -dz/2 -dz/2 0 dz/2 dz/2 0;
   dy/2 -dy/2 0 dy/2 -dy/2 0];
% thruster force [N]
F=1e-3;
% thruster force error (%) [-]
Te=0.03;
% thruster deactivation threshhold [N]
LL=1e-5;
% thruster activation threshhold [N]
UL=1e-4;

% REACTION WHEEL

% reaction wheel configuration matrix
ARW=[0 0; 1 0; 0 1];
% maximum torque of reaction wheels [Nm]
Mr_max=5e-3;
% reaction wheel inertia [kgm^2]
Irw=80e-6;

%--------------------------------------------------------------------------

Tp=[999.999926226686e-003    355.891478141766e-006   -144.526391847740e-006;
   -355.834344037688e-006    999.999858618000e-003    395.153008235765e-006;
    144.667003002503e-006   -395.101551630179e-006    999.999911483107e-003];

RRW=[R(:,3) R(:,6)];

pR=pinv(R);
pRRW=pinv(RRW);
pARW=pinv(ARW);

%--------------------------------------------------------------------------

dt_det=2450;
dt_po=850-100;

% Detumble

% initial S/C true anomaly [rad]
th=deg2rad(0); 
% S/C initial angular velocity [rad/s]
w0=[0.0754    0.0380    0.0579];
% initial S/C A_B/N orientation
A0=[1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
A0v=[A0(1,:) A0(2,:) A0(3,:)]';

th0=asin(A0(2,3));
th0=th0-pi*fix(th0/pi);

if abs(rad2deg(th0-pi/2))<10
    E=0;    % 313 system
else
    E=1;    % 312 system
end

kep=[a e i OM om th];   % kep=[a,e,i,OM,om,th] vector of keplerian elements

disp('Execute Detumble')

%%
close all
clc

t1=out.t;                   % time span
wx=out.w(:,1);              % x component of w [rad/s]
wy=out.w(:,2);              % y component of w [rad/s]
wz=out.w(:,3);              % z component of w [rad/s]

w1=[wx,wy,wz];              % S/C angual velocity [rad/s]
wn=out.wn;                  % angular velocity norm [rad/s]

A1=out.A;                   % A_B/N direction cosine matrix
M1=out.Ma;

ub=out.bu;                  % u body frame unit vector
vb=out.bv;                  % v body frame unit vector
wb=out.bw;                  % w body frame unit vector

R1=out.R;                    % S/C position [km]

figure
tiledlayout(3,1)
nexttile
plot(t1,M1(:,1),'linewidth',2)
ylabel('T_{c_x} [Nm]')
xlabel('Time [s]')
grid on
ylim([-5e-4 5e-4])
title('Control Torque T_x')
nexttile
plot(t1,M1(:,2),'linewidth',2)
ylabel('T_{c_y} [Nm]')
xlabel('Time [s]')
grid on
ylim([-7e-4 7e-4])
title('Control Torque T_y')
nexttile
plot(t1,M1(:,3),'linewidth',2)
ylabel('T_{c_z} [Nm]')
xlabel('Time [s]')
grid on
ylim([-8e-4 8e-4])
title('Control Torque T_z')

figure
plot(t1,wx,t1,wy,t1,wz,'linewidth',2)
grid on
title('Angular Velocity Components')
xlabel('Time [s]')
ylabel('\omega_i [rad/s]')
legend('\omega_x','\omega_y','\omega_z','location','best')

figure
plot(t1,wn,'linewidth',2)
xlabel('Time [s]')
ylabel('|\omega| [rad/s]')
grid on
title('Angular Velocity Norm')
legend('||\omega_x||','location','best')

hb_=[];
hn_=[];
hb=[];
hn=[];

for j=1:length(t1)
    hb_(:,j)=J*w1(j,:)';
    hn_(:,j)=A1(:,:,j)'*hb_(:,j);
    hb(j)=norm(hb_(:,j));
    hn(j)=norm(hn_(:,j));
end

figure
plot(t1,hb_(1,:),t1,hb_(2,:),t1,hb_(3,:),t1,hb,'--')
title('Angular Momentum, Body Frame')
legend('h_{B_x}','h_{B_y}','h_{B_z}','|h_B|','location','best')

figure
plot(t1,hn_(1,:),t1,hn_(2,:),t1,hn_(3,:),t1,hn,'--')
title('Angular Momentum, Inertial Frame')
legend('h_{N_x}','h_{N_y}','h_{N_z}','|h_N|','location','best')

figure
fig=gcf;
nplot=150;
hold on
plot3(R1(:,1),R1(:,2),R1(:,3),'k');
for j=1:fix(length(t1)/nplot)+1:length(t1)
    x=quiver3(R1(j,1),R1(j,2),R1(j,3),4*H*ub(:,1,j),4*H*ub(:,2,j),4*H*ub(:,3,j),'r');
    y=quiver3(R1(j,1),R1(j,2),R1(j,3),4*H*vb(:,1,j),4*H*vb(:,2,j),4*H*vb(:,3,j),'b');
    z=quiver3(R1(j,1),R1(j,2),R1(j,3),4*H*wb(:,1,j),4*H*wb(:,2,j),4*H*wb(:,3,j),'g');
    rr=plot3([0, R1(j,1)],[0, R1(j,2)],[0, R1(j,3)],'k');
    legend('orb','x','y','z')
    view([0.5 1 0.5])
    xlim([-1.2*a 1.2*a])
    ylim([-1.2*a 1.2*a])
    zlim([-1.2*a 1.2*a])
    xlabel('\gamma')
    ylabel('y')
    zlabel('z')
    drawnow
    delete(x)
    delete(y)
    delete(z)
    delete(rr)
    if ~ishghandle(fig)
        break
    end
end

figure
fig=gcf;
for j=1:fix(length(t1)/nplot)+1:length(t1)
    hold on
    x=quiver3(0,0,0,ub(:,1,j),ub(:,2,j),ub(:,3,j),'r');
    y=quiver3(0,0,0,vb(:,1,j),vb(:,2,j),vb(:,3,j),'b');
    z=quiver3(0,0,0,wb(:,1,j),wb(:,2,j),wb(:,3,j),'g');
    hold off
    view([0.5 1 0.5])
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    legend('x','y','z')
    drawnow
    delete(x)
    delete(y)
    delete(z)
    if ~ishghandle(fig)
        break
    end
end

% Final attitude parameters

tf=dt_det;
Af=A1(:,:,end);
thf=out.th(end);
wf=w1(end,:);

%% Slew

% Simulation initial time
t0=tf;
% S/C initial angular velocity [rad/s]
w0=wf;
% initial S/C A_B/N orientation
A0=Af;
A0v=[A0(1,:) A0(2,:) A0(3,:)]';

th0=asin(A0(2,3));
th0=th0-pi*fix(th0/pi);

if abs(rad2deg(th0-pi/2))<10
    E=0;    % 313 system
else
    E=1;    % 312 system
end

kep=[a e i OM om thf];   % kep=[a,e,i,OM,om,th] vector of keplerian elements

disp('Execute Slew')

%%
close all
clc

t2=out.t;                   % time span
wx=out.w(:,1);              % x component of w [rad/s]
wy=out.w(:,2);              % y component of w [rad/s]
wz=out.w(:,3);              % z component of w [rad/s]

w2=[wx,wy,wz];              % S/C angual velocity [rad/s]

A2=out.A;                   % A_B/N direction cosine matrix
M2=out.Ma;

ax=rad2deg(out.alphas(:,1));
ay=rad2deg(out.alphas(:,2));
az=rad2deg(out.alphas(:,3));

ub=out.bu;                  % u body frame unit vector
vb=out.bv;                  % v body frame unit vector
wb=out.bw;                  % w body frame unit vector

R2=out.R;                   % S/C position [km]

figure
tiledlayout(3,1)
nexttile
plot(t2,ax,'linewidth',2)
xlabel('Time [s]')
ylabel('\alpha_x [deg]')
grid on
title('Control Parameter \alpha_x')
nexttile
plot(t2,ay,'linewidth',2)
xlabel('Time [s]')
ylabel('\alpha_y [deg]')
grid on
title('Control Parameter \alpha_y')
nexttile
plot(t2,az,'linewidth',2)
xlabel('Time [s]')
ylabel('\alpha_z [deg]')
grid on
title('Control Parameter \alpha_z')

figure
tiledlayout(3,1)
nexttile
plot(t2,M2(:,1),'linewidth',2)
ylabel('T_{c_x} [Nm]')
xlabel('Time [s]')
grid on
ylim([-5e-4 5e-4])
title('Control Torque T_x')
nexttile
plot(t2,M2(:,2),'linewidth',2)
ylabel('T_{c_y} [Nm]')
xlabel('Time [s]')
grid on
ylim([-7e-4 7e-4])
title('Control Torque T_y')
nexttile
plot(t2,M2(:,3),'linewidth',2)
ylabel('T_{c_z} [Nm]')
xlabel('Time [s]')
grid on
ylim([-8e-4 8e-4])
title('Control Torque T_z')

figure
plot(t2,wx,t2,wy,t2,wz,'linewidth',2)
ylabel('\omega_i [rad/s]')
xlabel('Time [s]')
grid on
title('Angular Velocity Components')
legend('\omega_x','\omega_y','\omega_z','location','best')

hb_=[];
hn_=[];
hb=[];
hn=[];

for j=1:length(t2)
    hb_(:,j)=J*w2(j,:)';
    hn_(:,j)=A2(:,:,j)'*hb_(:,j);
    hb(j)=norm(hb_(:,j));
    hn(j)=norm(hn_(:,j));
    kin_(j)=0.5*w2(j,:)*hb_(:,j);
end

figure
plot(t2,hb_(1,:),t2,hb_(2,:),t2,hb_(3,:),t2,hb,'--')
title('Angular Momentum, Body Frame')
legend('h_{B_x}','h_{B_y}','h_{B_z}','|h_B|','location','best')

figure
plot(t2,hn_(1,:),t2,hn_(2,:),t2,hn_(3,:),t2,hn,'--')
title('Angular Momentum, Inertial Frame')
legend('h_{N_x}','h_{N_y}','h_{N_z}','|h_N|','location','best')

lat=rad2deg(out.lat_lambda);
lon=rad2deg(out.long_phi);

figure
title('Nadir Pointing')
hold on
plot(lat,lon,'-')
plot(lat(end),lon(end),'xr','linewidth',3)
xlabel('Lat [deg]')
ylabel('Lon [deg]')
xlim([-5 5])
ylim([-5 5])
hold off

figure
fig=gcf;
nplot=250;
hold on
plot3(R2(:,1),R2(:,2),R2(:,3),'k');
for j=1:fix(length(t2)/nplot)+1:length(t2)
    x=quiver3(R2(j,1),R2(j,2),R2(j,3),4*H*ub(:,1,j),4*H*ub(:,2,j),4*H*ub(:,3,j),'r');
    y=quiver3(R2(j,1),R2(j,2),R2(j,3),4*H*vb(:,1,j),4*H*vb(:,2,j),4*H*vb(:,3,j),'b');
    z=quiver3(R2(j,1),R2(j,2),R2(j,3),4*H*wb(:,1,j),4*H*wb(:,2,j),4*H*wb(:,3,j),'g');
    rr=plot3([0, R2(j,1)],[0, R2(j,2)],[0, R2(j,3)],'k');
    legend('orb','x','y','z')
    view([0.5 1 0.5])
    xlim([-1.2*a 1.2*a])
    ylim([-1.2*a 1.2*a])
    zlim([-1.2*a 1.2*a])
    xlabel('\gamma')
    ylabel('y')
    zlabel('z')
    drawnow
    delete(x)
    delete(y)
    delete(z)
    delete(rr)
    if ~ishghandle(fig)
        break
    end
end

figure
fig=gcf;
for j=1:fix(length(t2)/nplot)+1:length(t2)
    hold on
    x=quiver3(0,0,0,ub(:,1,j),ub(:,2,j),ub(:,3,j),'r');
    y=quiver3(0,0,0,vb(:,1,j),vb(:,2,j),vb(:,3,j),'b');
    z=quiver3(0,0,0,wb(:,1,j),wb(:,2,j),wb(:,3,j),'g');
    hold off
    view([0.5 1 0.5])
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    legend('x','y','z')
    drawnow
    delete(x)
    delete(y)
    delete(z)
    if ~ishghandle(fig)
        break
    end
end

% Final attitude parameters

tf=dt_po+t0;
Af=A2(:,:,end);
thf=out.th(end);
wf=w2(end,:);

%% Pointing

% Simulation initial time
t0=tf;
% initial S/C true anomaly [rad]
th=deg2rad(thf); 
% S/C initial angular velocity [rad/s]
w0=wf;
% initial S/C A_B/N orientation
A0=Af;
A0v=[A0(1,:) A0(2,:) A0(3,:)]';

th0=asin(A0(2,3));
th0=th0-pi*fix(th0/pi);

if abs(rad2deg(th0-pi/2))<10
    E=0;    % 313 system
else
    E=1;    % 312 system
end

kep=[a e i OM om thf];   % kep=[a,e,i,OM,om,th] vector of keplerian elements

%--------------------------------------------------------------------------

% initial conditions of linearized control

Aln0=Rz(thf)*Rx(i);             % local frame initial condition
Abl0=A0*Aln0';                  % attitude error initial condition
ax0=0.5*(Abl0(2,3)-Abl0(3,2));  % alpha_x initial condition [rad]
ay0=0.5*(Abl0(3,1)-Abl0(1,3));  % alpha_y initial condition [rad]  
az0=0.5*(Abl0(1,2)-Abl0(2,1));  % alpha_z initial condition [rad]
adx0=w0(1)-ay0*n;               % alpha_dot_x initial condition [rad/s]
ady0=w0(2)+ax0*n;               % alpha_dot_y initial condition [rad/s]
adz0=w0(3)-n;                   % alpha_dot_z initial condition [rad/s]
a0=[adx0 ady0 adz0 ax0 ay0 az0];

disp('Execute Pointing')

%% Idle pointing
close all
clc

t3=out.t;                   % time span
wx=out.w(:,1);              % x component of w [rad/s]
wy=out.w(:,2);              % y component of w [rad/s]
wz=out.w(:,3);              % z component of w [rad/s]

w3=[wx,wy,wz];              % S/C angual velocity [rad/s]

A3=out.A;                   % A_B/N direction cosine matrix
M3=out.Ma;

ax=rad2deg(out.alphas(:,1));
ay=rad2deg(out.alphas(:,2));
az=rad2deg(out.alphas(:,3));

ub=out.bu;                  % u body frame unit vector
vb=out.bv;                  % v body frame unit vector
wb=out.bw;                  % w body frame unit vector

R3=out.R;                   % S/C position [km]

figure
tiledlayout(3,1)
nexttile
plot(t3,ax,'linewidth',2)
xlabel('Time [s]')
ylabel('\alpha_x [deg]')
grid on
title('Control Parameter \alpha_x')
nexttile
plot(t3,ay,'linewidth',2)
xlabel('Time [s]')
ylabel('\alpha_y [deg]')
grid on
title('Control Parameter \alpha_y')
nexttile
plot(t3,az,'linewidth',2)
xlabel('Time [s]')
ylabel('\alpha_z [deg]')
grid on
title('Control Parameter \alpha_z')

figure
tiledlayout(3,1)
nexttile
plot(t3,M3(:,1),'linewidth',2)
ylabel('T_{c_x} [Nm]')
xlabel('Time [s]')
grid on
ylim([-5e-4 5e-4])
title('Control Torque T_x')
nexttile
plot(t3,M3(:,2),'linewidth',2)
ylabel('T_{c_y} [Nm]')
xlabel('Time [s]')
grid on
% ylim([-7e-4 7e-4])
title('Control Torque T_y')
nexttile
plot(t3,M3(:,3),'linewidth',2)
ylabel('T_{c_z} [Nm]')
xlabel('Time [s]')
grid on
% ylim([-8e-4 8e-4])
title('Control Torque T_z')

figure
plot(t3,wx,t3,wy,t3,wz,'linewidth',2)
ylabel('\omega_i [rad/s]')
xlabel('Time [s]')
title('Angular Velocity Components')
legend('\omega_x','\omega_y','\omega_z','location','best')

hb_=[];
hn_=[];
hb=[];
hn=[];

for j=1:length(t3)
    hb_(:,j)=J*w3(j,:)';
    hn_(:,j)=A3(:,:,j)'*hb_(:,j);
    hb(j)=norm(hb_(:,j));
    hn(j)=norm(hn_(:,j));
    kin_(j)=0.5*w3(j,:)*hb_(:,j);
end

figure
plot(t3,hb_(1,:),t3,hb_(2,:),t3,hb_(3,:),t3,hb,'--')
title('Angular Momentum, Body Frame')
legend('h_{B_x}','h_{B_y}','h_{B_z}','|h_B|','location','best')

figure
plot(t3,hn_(1,:),t3,hn_(2,:),t3,hn_(3,:),t3,hn,'--')
title('Angular Momentum, Inertial Frame')
legend('h_{N_x}','h_{N_y}','h_{N_z}','|h_N|','location','best')

lat=rad2deg(out.lat_lambda);
lon=rad2deg(out.long_phi);

figure
title('Nadir Pointing')
hold on
plot(lat,lon,'-')
plot(lat(end),lon(end),'xr','linewidth',3)
plot(lat(1),lon(1),'og','linewidth',3)
xlabel('\lambda [deg]')
ylabel('\phi [deg]')
xlim([-nad nad])
ylim([-nad nad])
hold off

figure
fig=gcf;
nplot=250;
hold on
plot3(R3(:,1),R3(:,2),R3(:,3),'k');
for j=1:fix(length(t3)/nplot)+1:length(t3)
    x=quiver3(R3(j,1),R3(j,2),R3(j,3),4*H*ub(:,1,j),4*H*ub(:,2,j),4*H*ub(:,3,j),'r');
    y=quiver3(R3(j,1),R3(j,2),R3(j,3),4*H*vb(:,1,j),4*H*vb(:,2,j),4*H*vb(:,3,j),'b');
    z=quiver3(R3(j,1),R3(j,2),R3(j,3),4*H*wb(:,1,j),4*H*wb(:,2,j),4*H*wb(:,3,j),'g');
    rr=plot3([0, R3(j,1)],[0, R3(j,2)],[0, R3(j,3)],'k');
    legend('orb','x','y','z')
    view([0.5 1 0.5])
    xlim([-1.2*a 1.2*a])
    ylim([-1.2*a 1.2*a])
    zlim([-1.2*a 1.2*a])
    xlabel('\gamma')
    ylabel('y')
    zlabel('z')
    drawnow
    delete(x)
    delete(y)
    delete(z)
    delete(rr)
    if ~ishghandle(fig)
        break
    end
end

figure
fig=gcf;
for j=1:fix(length(t3)/nplot)+1:length(t3)
    hold on
    x=quiver3(0,0,0,ub(:,1,j),ub(:,2,j),ub(:,3,j),'r');
    y=quiver3(0,0,0,vb(:,1,j),vb(:,2,j),vb(:,3,j),'b');
    z=quiver3(0,0,0,wb(:,1,j),wb(:,2,j),wb(:,3,j),'g');
    hold off
    view([0.5 1 0.5])
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    legend('x','y','z')
    drawnow
    delete(x)
    delete(y)
    delete(z)
    if ~ishghandle(fig)
        break
    end
end

%% Controlled Pointing
close all
clc

t4=out.t;                   % time span
wx=out.w(:,1);              % x component of w [rad/s]
wy=out.w(:,2);              % y component of w [rad/s]
wz=out.w(:,3);              % z component of w [rad/s]

w4=[wx,wy,wz];              % S/C angual velocity [rad/s]

A4=out.A;                   % A_B/N direction cosine matrix
M4=out.Ma;

ax=rad2deg(out.alphas(:,1));
ay=rad2deg(out.alphas(:,2));
az=rad2deg(out.alphas(:,3));

ub=out.bu;                  % u body frame unit vector
vb=out.bv;                  % v body frame unit vector
wb=out.bw;                  % w body frame unit vector

figure
tiledlayout(3,1)
nexttile
plot(t4,ax,'linewidth',2)
xlabel('Time [s]')
ylabel('\alpha_x [deg]')
grid on
title('Control Parameter \alpha_x')
nexttile
plot(t4,ay,'linewidth',2)
xlabel('Time [s]')
ylabel('\alpha_y [deg]')
grid on
title('Control Parameter \alpha_y')
nexttile
plot(t4,az,'linewidth',2)
xlabel('Time [s]')
ylabel('\alpha_z [deg]')
grid on
title('Control Parameter \alpha_z')

figure
tiledlayout(3,1)
nexttile
plot(t4,M4(:,1),'linewidth',2)
ylabel('T_{c_x} [Nm]')
xlabel('Time [s]')
grid on
ylim([-5e-4 5e-4])
title('Control Torque T_x')
nexttile
plot(t4,M4(:,2),'linewidth',2)
ylabel('T_{c_y} [Nm]')
xlabel('Time [s]')
grid on
ylim([-7e-4 7e-4])
title('Control Torque T_y')
nexttile
plot(t4,M4(:,3),'linewidth',2)
ylabel('T_{c_z} [Nm]')
xlabel('Time [s]')
grid on
ylim([-8e-4 8e-4])
title('Control Torque T_z')

figure
plot(t4,wx,t4,wy,t4,wz,'linewidth',2)
ylabel('\omega_i [rad/s]')
xlabel('Time [s]')
title('Angular Velocity Components')
legend('\omega_x','\omega_y','\omega_z','location','best')

lat=rad2deg(out.lat_lambda);
lon=rad2deg(out.long_phi);

figure
title('Nadir Pointing')
hold on
plot(lat,lon,'-')
plot(lat(end),lon(end),'xr','linewidth',3)
plot(lat(1),lon(1),'og','linewidth',3)
xlabel('\lambda [deg]')
ylabel('\phi [deg]')
xlim([-(nad+5) (5+nad)])
ylim([-(nad+13) (13+nad)])
hold off

%% End Plot

t=[t1' t2' t3'];
w=[w1; w2; w3;];

ub1=A1(1,:,:);
uv1=A1(2,:,:);
uw1=A1(3,:,:);

ub2=A2(1,:,:);
uv2=A2(2,:,:);
uw2=A2(3,:,:);

ub3=A3(1,:,:);
uv3=A3(2,:,:);
uw3=A3(3,:,:);

A=cat(3,A1,A2,A3);

ub=cat(3,ub1,ub2,ub3);
uv=cat(3,uv1,uv2,uv3);
uw=cat(3,uw1,uw2,uw3);

R_sc=[R1; R2; R3];

Ma=[M1; M2; M3];

figure
fig=gcf;
nplot=350;
hold on
plot3(R_sc(:,1),R_sc(:,2),R_sc(:,3),'k');
for j=1:fix(length(t)/nplot)+1:length(t)
    x=quiver3(R_sc(j,1),R_sc(j,2),R_sc(j,3),4*H*ub(:,1,j),4*H*ub(:,2,j),4*H*ub(:,3,j),'r');
    y=quiver3(R_sc(j,1),R_sc(j,2),R_sc(j,3),4*H*uv(:,1,j),4*H*uv(:,2,j),4*H*uv(:,3,j),'b');
    z=quiver3(R_sc(j,1),R_sc(j,2),R_sc(j,3),4*H*uw(:,1,j),4*H*uw(:,2,j),4*H*uw(:,3,j),'g');
    rr=plot3([0, R_sc(j,1)],[0, R_sc(j,2)],[0, R_sc(j,3)],'k');
    legend('orb','x','y','z')
    view([0.5 1 0.5])
    xlim([-1.2*a 1.2*a])
    ylim([-1.2*a 1.2*a])
    zlim([-1.2*a 1.2*a])
    xlabel('\gamma')
    ylabel('y')
    zlabel('z')
    drawnow
    delete(x)
    delete(y)
    delete(z)
    delete(rr)
    if ~ishghandle(fig)
        break
    end
end

figure
hold on
plot(t,w(:,1),t,w(:,2),t,w(:,3),'linewidth',2)
plot([t1(end) t1(end)],[-0.1 0.1],'r--')
plot([t2(end) t2(end)],[-0.1 0.1],'r--')
hold off
xlabel('Time [s]')
legend('\omega_x','\omega_y','\omega_z')
title('Angular Velocity \omega')

for j=1:length(t)
    a11(j)=A(1,1,j);
    a12(j)=A(1,2,j);
    a13(j)=A(1,3,j);
    a21(j)=A(2,1,j);
    a22(j)=A(2,2,j);
    a23(j)=A(2,3,j);
    a31(j)=A(3,1,j);
    a32(j)=A(3,2,j);
    a33(j)=A(3,3,j);
end

figure
plot(t,a11,t,a12,t,a13,t,a21,t,a22,t,a23,t,a31,t,a32,t,a33)
title('Attitude DCM')

figure
tiledlayout(3,1)
nexttile
hold on
plot(t,Ma(:,1),'linewidth',2)
plot([t1(end) t1(end)],[-5e-4 5e-4],'r--')
plot([t2(end) t2(end)],[-5e-4 5e-4],'r--')
hold off
ylabel('T_{c_x} [Nm]')
xlabel('Time [s]')
grid on
ylim([-5e-4 5e-4])
title('Control Torque T_x')
nexttile
hold on
plot(t,Ma(:,2),'linewidth',2)
plot([t1(end) t1(end)],[-7e-4 7e-4],'r--')
plot([t2(end) t2(end)],[-7e-4 7e-4],'r--')
hold off
ylabel('T_{c_y} [Nm]')
xlabel('Time [s]')
grid on
ylim([-7e-4 7e-4])
title('Control Torque T_y')
nexttile
hold on
plot(t,Ma(:,3),'linewidth',2)
plot([t1(end) t1(end)],[-8e-4 8e-4],'r--')
plot([t2(end) t2(end)],[-8e-4 8e-4],'r--')
hold off
ylabel('T_{c_z} [Nm]')
xlabel('Time [s]')
grid on
ylim([-8e-4 8e-4])
title('Control Torque T_z')
