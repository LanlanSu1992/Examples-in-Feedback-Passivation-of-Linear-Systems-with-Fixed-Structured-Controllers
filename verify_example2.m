clear all;
syms s;

% kp= 0.516
% ki=0.669
%kd=-1;
 kp= 1
 ki=1

D=(s-1)*(s-2);
double(solve(D==0))
N=(s+2)*(s+3);



C=kp+ki/(s+1)%+kd*s;

G=(N/D)/(1+(N/D)*C);
[n,d] = numden(simplify(G));

[coeffd,s1]=coeffs(d, s)
[coeffn,s2]=coeffs(n, s)

[A,B,C,D]=tf2ss([ coeffn ],coeffd);
A=double(A);
B=double(B);
C=double(C);
D=double(D);
n=max(size(A));
m=1

P=sdpvar(n);
sdpvar xi;

M=[A.'*P+P*A-xi*C.'*C P*B-0.5*C.'-C.'*xi*D;(P*B-0.5*C.'-C.'*xi*D).' -D.'*xi*D-D];
CON=[P>=0,M<=0];
obj=xi;
sol=optimize(CON,obj);
xi_value=value(xi)

% P=sdpvar(n);
% sdpvar nu;
% S=0.5.*C.'*eye(m);
% R=D.'*0.5+0.5*D-nu;
% M=[A.'*P+P*A P*B-S;B.'*P.'-S.' -R];
% CON=[P>=0,M<=0];
% obj=-nu;
% sol=optimize(CON,obj);
% nu_value=value(nu)

