clc, clear all, close all
addpath("data\")
addpath("./functions")

alpha = load('data/alpha_a3_f1.mat').alpha(:,2);
theta = load('data/theta_a3_f1.mat').theta(:,2);
% load("theta_a1_f1.mat");
t = (0:0.01:40);
u1 = 0.01*sin(t);
alpha_dot = [0;100*(alpha(2:end)-alpha(1:end-1))];
theta_dot = [0;100*(theta(2:end)-theta(1:end-1))];
%%


[Abar,Bbar,C,D,x0, J, H] = pem(theta,A0,B0,C0,D0,x00,y,u,lambda,maxiter)


%%



function [Abar,Bbar,C,D,x0] = theta2matrices(theta,Asize,Bsize,Csize,Dsize,xsize)

r_Abar = Asize(1);
c_Abar = Asize(2);
r_Bbar = Bsize(1);
c_Bbar = Bsize(2);
r_C = Csize(1);
c_C = Csize(2);
r_D = Dsize(1);
c_D = Dsize(2);
l_x0 = xsize;

for i = 1:r_Abar
    for j = 1:c_Abar
        Abar(i,j) = theta(i+(j-1)*r_Abar);
    end
end
for i = 1:r_Bbar
    for j = 1:c_Bbar
        Bbar(i,j) = theta(i+(j-1)*r_Bbar+r_Abar*r_Bbar);
    end
end

for i = 1:r_C
    for j = 1:c_C
        C(i,j) = theta(i+(j-1)*r_C+r_Abar*c_Abar+r_Bbar*c_Bbar);
    end
end

for i = 1:r_D
    for j = 1:c_D
        D(i,j) = theta(i+(j-1)*r_D+r_Abar*c_Abar+r_Bbar*c_Bbar +r_C*c_C);
    end
end
for i = 1:l_x0
    x0(i) = theta(i +r_Abar*c_Abar+r_Bbar*c_Bbar +r_C*c_C +r_D*c_D);
end
x0 = x0';

end