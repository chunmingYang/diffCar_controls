function J = diffCar_MPC_costFunction(X,U,Xref,ts,Q,R)
    J = 0;
    for i = 1 : length(U)
%         J = J + transpose(Xref - X)*Q*(Xref - X) + u*R*u; % since in this
%         case we don't need to optimize inputs/controls so we don't have
%         "u*R*u" this term, and this case Xref is fixed during the "i = 1
%         : length(U)" receding horizon period
        J = J + transpose(Xref - X)*Q*(Xref - X); 

        % iteration to update states
        A = [0 0 0;
             0 0 0;
             0 0 0];
        B = [cos(X(3)) -sin(X(3)) 0;
             sin(X(3))  cos(X(3)) 0;
             0          0         1];
        Xdot = A*X + B*U(:,i);
        X = X + ts*Xdot;

    end
end