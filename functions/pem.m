function [Abar,Bbar,C,D,x0, J, H] = pem(theta,A0,B0,C0,D0,x00,y,u,lambda,maxiter)
% Instructions:
% Implement your Prediction Error Method for the Output Error system here.
% Use the following function inputs and outputs.
%
% Function INPUT
% theta  Paramter vector (size: depends on your mapping choice)
% A0 Initial guess for system matrix A (matrix of size n x n)
% B0 Initial guess for system matrix B (matrix of size n x m)
% C0 Initial guess for system matrix C (matrix of size l x n)
% D0 Initial guess for system matrix D (matrix of size l x m)
% x00 Initial guess for initial state (vector of size n x one)
% u System input (matrix of size N x m)
% y System output (matrix of size N x l)
% lambda regularization parameter (scalar)
% maxiter Maximum number of iterations (scalar)
%
%
% Function OUTPUT
% Abar Estimate of system matrix A (matrix of size n x n)
% Bbar Estimate of system matrix B (matrix of size n x m)
% C Estimate of system matrix C (matrix of size l x n)
% D Estimate of system matrix D (matrix of size l x m)
% x0 Estimate of initial state (vector of size n x one)



% YOUR CODE HERE
% size parameter extraction
[n,~] = size(A0);
[~,m] = size(B0);
[l,~] = size(C0);
[N,~] = size(u);
p= length(theta);
% initial estimates
Abar = A0;
Bbar = B0;
C = C0;
D = D0;
x_hat(:,1)=x00;
% partial derivative calculation for theta_i
for i = 1:p
    Abar_dot(:,:,i) = zeros([n,n]);
    Bbar_dot(:,:,i)=zeros([n,1]);
    x0_dot(:,:,i) = zeros([n,1]);
    if i<6
        Abar_dot(i,1,i)=1;
    elseif i<9
        Bbar_dot(i-5,1,i) = 1;
    else
        x0_dot(i-8,1,i) = 1;
    end
end
C_dot = 0;
D_dot = 0;
K_dot = 0;

% Gauss-Newton

% initialize
converged = false;
maxiter_reached = false;
iter = 1;
theta_new = theta;

while ~converged && ~maxiter_reached

    % YOUR CODE HERE
    x_hat_dot = zeros(length(theta),n,length(u));
    for i=1:p
        for j=1:N-1
            x_hat(:,j+1) = Abar*x_hat(:,j)+Bbar*u(j,:);
            y_hat(j) = C*x_hat(:,j);
            x_hat_dot(i,:,j+1)=Abar*x_hat_dot(i,:,j)'+Abar_dot(:,:,i)*x_hat(:,j)+Bbar_dot(:,:,i).*u(j,:);
            y_hat_dot(i,j) = C*x_hat_dot(i,:,j)';
            E_n(j)=(y(j)-y_hat(j));
        end
        psi(i,:) = -y_hat_dot(i,:);
        J(i) = 2/N*psi(i,:)*E_n';
    end
    H = 2/N*(psi*psi');
    theta = theta-pinv(H+lambda*eye(size(H)))*J';
    [Abar,Bbar,C,D,x0] = theta2matrices(theta);
    % Check convergence
    if norm(theta_new - theta) < 1e-3 
        converged = true; 
    elseif iter == maxiter 
        maxiter_reached = true; 
        warning('Maximum iterations reached'); 
    end
    theta = theta_new; 
    iter = iter + 1;
end

     
end
