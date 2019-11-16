clear all;
syms z;
D=(z-2);
double(solve(D==0))
N=z;

syms k_p k_i 'real'
C=k_p+k_i/(z-0.5);





c1=1-k_p;
c2=k_p-0.1;
c3=2-k_i;
c4=k_i-1;


G=(N/D)/(1+(N/D)*C);
[n,d] = numden(simplify(G));
[coeffd,z1]=coeffs(d, z)
[coeffn,z2]=coeffs(n, z)
syms y 'real'
N_c=subs(n,z,(1-y^2+sqrt(-1)*2*y)/(1+y^2));
D_c=subs(d,z,(1-y^2+sqrt(-1)*2*y)/(1+y^2));


P_ne=simplify(real(N_c*(1+y^2)^2));
P_no=simplify(imag(N_c*(1+y^2)^2));

P_de=simplify(real(D_c*(1+y^2)^2));
P_do=simplify(imag(D_c*(1+y^2)^2));

g=P_ne*P_de+P_no*P_do

% syms eps
% prog=sosprogram(y);
% prog=sosdecvar(prog,[eps k_p k_i]);
%      
%      
%       prog=sosineq(prog,c1);
%       prog=sosineq(prog,c2);
%       prog=sosineq(prog,c3);
%       prog=sosineq(prog,c4);
% 
%       
%      prog=sosineq(prog,g-eps);
%      prog=sossetobj(prog,-eps);
%      prog=sossolve(prog);
%      
%      eps=sosgetsol(prog,eps)
%      k_p1=sosgetsol(prog,k_p)
%      k_i1=sosgetsol(prog,k_i)

%      

syms nu e;
nu_upper=1;
nu_lower=0;
error=1;

while(error>1e-4)
    nu=1/2*(nu_upper+nu_lower);
    prog1=sosprogram(y);
    prog1=sosdecvar(prog1,[e k_p k_i]);
    M=[g sqrt(nu)*P_de sqrt(nu)*P_do;
    sqrt(nu)*P_de 1 0;
    sqrt(nu)*P_do 0 1];
    prog1=sosineq(prog1,c1);
    prog1=sosineq(prog1,c2);
    prog1=sosineq(prog1,c3);
    prog1=sosineq(prog1,c4);
     prog1=sosmatrixineq(prog1,M-e*eye(3));
     prog1=sossetobj(prog1,-e);
     prog1=sossolve(prog1);
     e_v=sosgetsol(prog1,e)
     if e_v>0

         nu_lower=nu;
     else
         nu_upper=nu;
     end
     error=nu_upper-nu_lower;
end
nu


     k_p1=double(sosgetsol(prog1,k_p))
     k_i1=double(sosgetsol(prog1,k_i))