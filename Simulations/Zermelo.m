clear all
clc

% Time Variables
syms t tf_o %s

% Parameters
V = 600 %m/s   Aircraft Speed
w = 60  %m/s   Wind Speed

po = -1;    % Lambdazero

theta_trim = 0.046492; % Trimming Angle

% Variables definitions
syms x1(t) x2(t) p1(t) p2(t) u(t)

% State equations
Dx1 = V*cos(u(t))
Dx2 = V*sin(u(t)) + w

% Cost function (inside integral)
syms J(t);
J = 1

% Boundary conditions
cond_1 = x1(0) == 0;
cond_2 = x2(0) == 0;
cond_3 = x1(tf_o) == 1;
cond_4 = x2(tf_o) == 0;

% Hamiltonian
syms H(t)
H = po*J + p1*Dx1 + p2*Dx2

% Costate equations
cost_1 = diff(p1,t) == -diff(H,x1);
cost_2 = diff(p2,t) == -diff(H,x2);

% Control equation
control_eqn = diff(H,u) == 0

% Solve costate
[C1,C2] = dsolve(cost_1,cost_2)

% Substitue solution of costate equations in control equation
assume(sin(u(t)) ~= 0);
assume(cos(u(t)) ~= 0);
assume(tan(u(t)) ~= 0);

control_eqn = subs(control_eqn/V,[p1,p2],[C1,C2]);

% Solve control equation
sol_step1 = isolate(control_eqn, C1);
sol_step2 = subs(sol_step1,sin(u(t))/cos(u(t)),tan(u(t)));
sol = isolate(sol_step2,u(t));

% Dropping time dependence of the u(t)
syms u_t
sol = subs(sol,u(t),u_t)
Dx1 = subs(Dx1,u(t),u_t);
Dx2 = subs(Dx2,u(t),u_t);

% Solutions
SOL_A1 = x1 == int(Dx1,t)
SOL_A2 = x2 == int(Dx2,t)

% Solve with Boundary Conditions
B_A1 = subs(SOL_A1,[x1,t],[rhs(cond_3),tf_o]);
B_A2 = subs(SOL_A2,x2,rhs(cond_4))/t;

% Find optimal control u_opt
u_opt = rhs(isolate(B_A2,u_t));

% Find optimal cos(u) and sin(u)
sin_opt = isolate(B_A2,sin(u_t));
cos_opt = simplify(isolate(rewrite(sin_opt^2,"cos"),cos(u_t)));

% Substitute to find optimal final time expression tf*
tf_opt1 = solve(B_A1,tf_o);
tf_opt = subs(tf_opt1,cos(u_t),rhs(cos_opt));

% Analytical solutions
u_opt = double(u_opt)
tf_o = double(tf_opt)

%-------------------------------------------------------------------
% Plots
% Optimal state variables evolution
syms x1_opt(t) x2_opt(t)
x1_opt = V*cos(u_opt)*t
x2_opt = (V*sin(u_opt)+ w)*t

hold on
x = linspace(1,10,200);

% Position along x1
figure(1)
fplot(x1_opt,[0 5],'b','LineWidth',2)
set(gca,'FontSize',14)
title('Longitudinal Position (m)','FontSize', 18)
xlabel('Time (s)','FontSize', 10) 
ylabel('x1(t)','FontSize', 10) 

% Position along x2
figure(2)
fplot(x2_opt + 15000,[0 5],'b','LineWidth',2)
set(gca,'FontSize',14)
title('Vertical Position (m)','FontSize', 18)
xlabel('Time (s)','FontSize', 10) 
ylabel('x2(t)','FontSize', 10) 
