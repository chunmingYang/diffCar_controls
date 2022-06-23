clc
clear all
close all

%% code description
% using LQR to do the trajectory tracking with diff car model that contains
% Vx, Vy, Omega as the inputs

%% model description
% Xdot = AX + BU
% Y = CX + DU
% states: X = [x; y; theta]
% inputs/controls: U = [vx; vy; omega]
% outputs: Y = [x; y; theta]
% reference: Yr = [xref; yref; thetaref]

% states/ref/controls initialization
X = [0;0;0];
Xdot = [0;0;0];
U = [0;0;0];
Yr = [0;0;0];

% dynamics initialization
A = [0 0 0;
    0 0 0;
    0 0 0];
B = [cos(X(3)) -sin(X(3)) 0;
    sin(X(3)) cos(X(3)) 0;
    0 0 1];
C = 1;
D = 0;

%% LQR parameters
Q = diag([1e3 1e3 1e0]);     % this can be tuned, according to X'*Q*X
R = 1;                       % this can be tuned, according to U'*R*U
ts = 0.001;                  % states update resolution
TS = 0.01;                   % step response resolution --- reference resolution
iter = TS/ts;                % iteration times within one single step response
ITER = 1000;                 % the whole simulation time = ITER*TS

%% reference trajectory generation
t0 = 0;
tf = TS*ITER;
y0 = 0;                      % position start position
yf = 0.5;                    % position end position

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
count = 1;

%% main loop
for I = 1 : ITER
    xref = a0 + a1*(I*TS) + a2*(I*TS)*(I*TS) + a3*(I*TS)*(I*TS)*(I*TS);
    yref = a0 + a1*(I*TS) + a2*(I*TS)*(I*TS) + a3*(I*TS)*(I*TS)*(I*TS);
    thetaref = 0;
    Yr = [xref; yref; thetaref];
    for i = 1 : iter

        % states update
        A = [0 0 0;
            0 0 0;
            0 0 0];
        B = [cos(X(3)) -sin(X(3)) 0;
            sin(X(3)) cos(X(3)) 0;
            0 0 1];
        C = 1;
        D = 0;
        Xdot = A*X + B*U;
        X = X + ts*Xdot;

        % optimization core
        P=are(A, B*inv(R)*B', C'*Q*C);
        g=inv(P*B*inv(R)*B'-A')*C'*Q*Yr;
        U = -inv(R)*B'*(P*X-g);

        % data saving
        xhis(count) = X(1);
        yhis(count) = X(2);
        thetahis(count) = X(3);
        xrefhis(count) = Yr(1);
        yrefhis(count) = Yr(2);
        thetarefhis(count) = Yr(3);
        count = count + 1;
        
    end
end

%% plots interpretation
t = [];
for i = 1:length(xhis)
    t(i) = i*ts;
end
subplot(3,1,1)
plot(t, xhis); hold on;
plot(t, xrefhis)
xlabel("x position")
legend("actual", "ref")

subplot(3,1,2)
plot(t, yhis); hold on;
plot(t, yrefhis)
xlabel("y position")
legend("actual", "ref")

subplot(3,1,3)
plot(t, thetahis); hold on;
plot(t, thetarefhis)
xlabel("yaw angle")
legend("actual", "ref")