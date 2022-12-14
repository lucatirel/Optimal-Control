%% SYSTEM DEFINITION
A = [0 0 1; 
    -1.1714*10^-13 -0.76788 0.9396;
    0 -2.251 -1.0462];
B = [0;
    -0.0016222;
    -0.16416];
C = [1 0 0];
D = 0;

x0 = [0.046492; 0.046492; 0]; %initial conditions

%% SYSTEM ANALYSIS
% System definition
sys = ss(A,B,C,D)
H = tf(sys)
SIZE = size(sys);
OPEN_LOOP = allmargin(H)

% Zeroes and Poles
P = pole(H)
Z = zero(H)

% Check reachability and observability
CTRB = ctrb(A,B);
rk_CTRB = rank(CTRB)
OBSV = obsv(A,C);
rk_OBSV = rank(OBSV)

%------------------------------------------------------
% METHOD 1: POLE PLACEMENT
% Desired closed-loop poles
PO = [-2+0.4j -2-0.4j -6]

% Pole placement gain
K1 = place(A,B,PO)

% Precompensator gain
G1 = -inv(C*inv(A-B*K1)*B)
%G1 = -33.9859 (wrong)

% Definition of closed loop system (Pole Placement)
sys_cl_1 = ss(A-B*K1,B*G1,C,D);
INFO_CLSYS_1 = stepinfo(sys_cl_1)
%------------------------------------------------------
% METHOD 2: LQR
% COST INDEX A
% Parameters definition
Q1 = [90 0 0;
      0 0 0;
      0 0 0]
R1 = 0.0001

% Optimal gain
K2 = lqr(A,B,Q1,R1)

% Find precompensator gain
G2 = -inv(C*inv(A-B*K2)*B)

% Definition of closed loop system
sys_cl_2 = ss(A-B*K2,B*G2,C,D);
INFO_CLSYS_2 = stepinfo(sys_cl_2)

%-------------------------------------------------------
% COST INDEX B
% Parameters definition
Q2 = [90 0 0;
    0 0 0;
    0 0 15]
R2 = 0.0001

% Optimal gain
K3 = lqr(A,B,Q2,R2)

% Precompensator gain
G3 = -inv(C*inv(A-B*K3)*B)

% Definition of closed loop system
sys_cl_3 = ss(A-B*K3,B*G3,C,D);
INFO_CLSYS_3 = stepinfo(sys_cl_3)

%---------------------------------------------------------
% METHOD 4: LQR with Minimal Control
% Parameters definition
% State Variable Weights
Q_u = [90 0 0;
       0 0 0;
       0 0 0]
% R values
r1 = 5;
r2 = 10;
r3 = 15;

% Gains
K_f = lqr(A,B,Q_u,r1);
K_s = lqr(A,B,Q_u,r2);
K_t = lqr(A,B,Q_u,r3);

% Precompensators 
G_f = -inv(C*inv(A-B*K_f)*B);
G_s = -inv(C*inv(A-B*K_s)*B);
G_t = -inv(C*inv(A-B*K_t)*B);

%PLOTS IN SIMULINK !
%---------------------------------------------------------
% METHOD 4: LQR with Bounded Control
% COST INDEX C
% Parameters definition
Q3 = [900 0 0;
      0 0 0;
      0 0 9]
R3 = 0.1

% Optimal gain
K4 = lqr(A,B,Q3,R3)
G4 = -inv(C*inv(A-B*K4)*B)

%PLOTS IN SIMULINK !
%---------------------------------------------------------
%% Plots
% Options 
opt = stepDataOptions("StepAmplitude",u_opt); % STEP AMPLITUDE = THETA OPTIMAL
bode_opt = bodeoptions('cstprefs');
bode_opt.MagVisible = 'off';

% Open Loop Analysis
figure(3)
subplot(2,2,1)
rlocus(sys) % Root Locus
subplot(2,2,2)
bodemag(sys) % Bode Magnitude
subplot(2,2,3)
bodeplot(sys,bode_opt) % Bode Phase
subplot(2,2,4)
step(sys,opt) % Step Response
ytickformat('%,.0f')

% Closed Loop with Pole Placement
figure(4)
[yl_1,t_1,xl_1] = initial(sys_cl_1,x0);
[yf_1,t_1,xf_1] = step(sys_cl_1,opt);
a_1 = yl_1(1:300);
b_1 = yf_1(1:300);
c_1 = a_1 + b_1;
plot(t_1,c_1)

title('Tracking the Optimal Pitch Angle using Pole Placement Method','FontSize', 16)
set(gca,'FontSize',16)
xlabel('Time','FontSize', 16) 
ylabel('Angle of Attack (Radians)','FontSize', 16) 

% Closed Loop with LQR and Q1,R1
figure(5)
[yl_2,t_2,xl_2] = initial(sys_cl_2,x0);
[yf_2,t_2,xf_2] = step(sys_cl_2,opt);
a_2 = yl_2;
b_2 = yf_2;
c_2 = a_2 + b_2;
plot(t_2,c_2)

title('Tracking the Optimal Pitch Angle using LQR with Q1 & R1','FontSize', 16)
set(gca,'FontSize',16)
xlabel('Time','FontSize', 16) 
ylabel('Angle of Attack (Radians)','FontSize', 16) 

% Closed Loop with LQR and Q2,R2
figure(6)
[yl_3,t_3,xl_3] = initial(sys_cl_3,x0);
[yf_3,t_3,xf_3] = step(sys_cl_3,opt);
a_3 = yl_3;
b_3 = yf_3;
c_3 = a_3 + b_3;
plot(t_3,c_3)

title('Tracking the Optimal Pitch Angle using LQR with Q2 & R2','FontSize', 16)
set(gca,'FontSize',16)
xlabel('Time','FontSize', 16) 
ylabel('Angle of Attack (Radians)','FontSize', 16) 

%--------------------------------------------------------------------------
% Notes
% Alternative Tuning Method
p = 50;
Q4 = p*C'*C;

%RiseTime — Time it takes for the response to rise from 10% to 90% of the steady-state response.
%SettlingTime — Time it takes for the error |y(t) - yfinal| between the response y(t) and the steady-state response yfinal to fall to within 2% of yfinal.
%Overshoot — Percentage overshoot, relative to yfinal).
%PeakTime — Time at which the peak value occurs.

