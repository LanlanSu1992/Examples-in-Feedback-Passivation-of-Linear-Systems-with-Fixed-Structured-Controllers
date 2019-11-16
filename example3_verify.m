clear all;
syms z;


D=z-2;
double(solve(D==0))
N=z;
kp= 0.1
ki=1.49
%kd=-1;
C=kp+ki/(z-0.5);%+k_d*s;



G=(N/D)/(1+(N/D)*C);




[n,d] = numden(simplify(G));

[coeffd,s1]=coeffs(d, z)
[coeffn,s2]=coeffs(n, z)

[A,B,C,D]=tf2ss([ coeffn 0],coeffd);
A=double(A);
B=double(B);
C=double(C);
D=double(D);
n=max(size(A));
m=1

P=sdpvar(n);
sdpvar nu;
S=0.5.*C.'*eye(m);
R=D.'*0.5+0.5*D-nu;
M=[A.'*P*A-P A.'*P*B-S;B.'*P.'*A-S.' -R+B.'*P*B];
CON=[P>=0,M<=0];
obj=-nu;
sol=optimize(CON,obj);
nu_value=value(nu)
