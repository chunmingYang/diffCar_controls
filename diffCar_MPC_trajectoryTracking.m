clc 
clear all
close all

%% MPC parameters
ts = 0.01;
T = 5;                % receding horizon can be tuned
iter = 200;           % iteration times, thus we have total_time = 200*0.01
n = 3;                % number of states (x, y, theta)
m = 3;                % number of inputs/controls (vx, vy, omega)
Q = diag([10 10 10]); % this can be tuned
R = zeros(m);         % inputs cost matrix since we don't have the target input trajectory to optimize then we set zero

%% initial guess for integrate problem

% states initialization
x = 0.1;
y = 0;
theta = 0;
X = [x; y; theta];

% states reference trajectory initialization
xref = 0;
yref = 0;
thetaref = 0;
Xref = [xref; yref; thetaref];

% inputs/controls initialization
U = zeros(m, T);

% diff car dynamics initialization
A = [0 0 0;
     0 0 0;
     0 0 0];
B = [cos(X(3)) -sin(X(3)) 0;
     sin(X(3))  cos(X(3)) 0;
     0          0         1];

%% cubic trajectory generation
t0 = 0;
tf = ts*iter;
y0 = x;    % trajectory start position
yf = 1.5;  % trajectory end position

a0 = (yf*t0*t0*(3*tf-t0) + y0*tf*tf*(tf-3*t0))/((tf-t0)*(tf-t0)*(tf-t0));
a1 = 6*t0*tf*(y0-yf)/((tf-t0)*(tf-t0)*(tf-t0));
a2 = 3*(t0+tf)*(yf-y0)/((tf-t0)*(tf-t0)*(tf-t0));
a3 = 2*(y0-yf)/((tf-t0)*(tf-t0)*(tf-t0));

%% data saving
xhis = [];
yhis = [];
thetahis = [];
xrefhis = [];
yrefhis = [];
thetarefhis = [];

%% MPC main loop
for i = 1 : iter
    Aieq = [];            % since no ineuqality constraint
    bieq = [];            % since no inequality constraint
    Aeq = [];             % since no equality constraint
    beq = [];             % since no equality constraint
    lb = -50*ones(m, T);  % inputs/controls low boundary
    ub = 50*ones(m, T);   % inputs/controls up boundary

    % cubic trajectory
    xref = a0 + a1*(i*ts) + a2*(i*ts)*(i*ts) + a3*(i*ts)*(i*ts)*(i*ts);
    yref = a0 + a1*(i*ts) + a2*(i*ts)*(i*ts) + a3*(i*ts)*(i*ts)*(i*ts);
    thetaref = 0;
    Xref = [xref; yref; thetaref];
    xrefhis(i) = xref;    % start from t=0.01
    yrefhis(i) = yref;
    thetarefhis(i) = thetaref;

    % optimization loop
    u = fmincon(@(U)diffCar_MPC_costFunction(X,U,Xref,ts,Q,R),U,Aieq,bieq,Aeq,beq,lb,ub);

    % iteration to update states
    A = [0 0 0;
         0 0 0;
         0 0 0];
    B = [cos(X(3)) -sin(X(3)) 0; ...
         sin(X(3))  cos(X(3)) 0;
         0          0         1];
    Xdot = A*X + B*u(:,1);  % core of MPC: using the first optimized "u"
    X = X + ts*Xdot;        % core of MPC: states iteration update

    % warm start
    U = repelem(u(:,1),1,T);
    
    % data saving
    xhis(i) = X(1);         % start from t=0.01
    yhis(i) = X(2);
    thetahis(i) = X(3)*180/3.14;
end

%% plot interpretation
t = [];
for i = 1:length(xhis)
    t(i) = i*ts;
end

subplot(3,1,1)
plot(t, xhis); hold on;
plot(t, xrefhis)
xlabel("x position")
legend("actual", "reference")

subplot(3,1,2)
plot(t, yhis); hold on;
plot(t, yrefhis)
xlabel("y position")
legend("actual", "reference")

subplot(3,1,3)
plot(t, thetahis); hold on;
plot(t, thetarefhis)
xlabel("yaw angle")
legend("actual", "reference")