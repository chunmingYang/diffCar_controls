clc
clear all
close all

%% code description
% LQR step response control using diff car model that contains Vx, Vy,
% Omega as the inputs inspried from (https://www.cnblogs.com/zhjblogs/p/13834109.html)

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
Yr = [0.5;0.5;0];

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
Q = diag([1e3 1e3 1]);    % this can be tuned, according to X'*Q*X
R = 1;                    % this can be tuned, according to U'*R*U
ts = 0.01;                % step time
iter = 500;               % iteration time

%% data saving
xhis = [];
yhis = [];
thetahis = [];

%% main loop
for i = 1 : iter

    % states update
    A = [0 0 0;
         0 0 0;
         0 0 0];
    B = [cos(X(3)) -sin(X(3)) 0;
         sin(X(3))  cos(X(3)) 0;
         0          0         1];
    C = 1;
    D = 0;
    Xdot = A*X + B*U;
    X = X + ts*Xdot;
    
    % optimization core
    P=are(A, B*inv(R)*B', C'*Q*C);
    g=inv(P*B*inv(R)*B'-A')*C'*Q*Yr;
    U = -inv(R)*B'*(P*X-g);
    
    % data saving
    xhis(i) = X(1);
    yhis(i) = X(2);
    thetahis(i) = X(3);
    
end

%% plots interpretation
t = [];
for i = 1:length(xhis)
    t(i) = i*ts;
end
subplot(3,1,1)
plot(t, xhis); hold on;
plot(t, Yr(1))
legend("actual", "ref")

subplot(3,1,2)
plot(t, yhis); hold on;
plot(t, Yr(2))
legend("actual", "ref")

subplot(3,1,3)
plot(t, thetahis); hold on;
plot(t, Yr(3))
legend("actual", "ref")






