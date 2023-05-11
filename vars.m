mr = 0.095;
Lr = 0.085;
mp = 0.024;
Lp = 0.024;
mh = 0.016;
Jh = 0.6*10^(-6);
Jm = 4*10^(-6);

Jstar = Jr*Jp-.25*mp^2*Lp^2*Lr*2;
A = 1/Jstar*[-Jp*br,.5*mp*Lp*Lr*bp,.25*mp^2*Lr*g;.5*mp*Lp*Lr*br,-Jr*bp,-.5*mp*g*Lp*Jr;0,Jstar,0];
B = 1/Jstar*[Jp;-.5*mp*Lp*Lr;0;0];
C = eye(3);
D=0;
