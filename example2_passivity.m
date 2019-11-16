clear all;% CT CASE STABULITY ANALYSS Part I
format long
syms s;
D=(s-1)*(s-2);

N=(s+2)*(s+3);

syms k_p k_i 'real'


C=k_p+k_i/(s+1);


%the range of PID

c1=1-k_p;
c2=k_p+0;
c3=1-k_i;
c4=k_i+0;

G=(N/D)/(1+(N/D)*C);
[n,d] = numden(simplify(G));

syms w 'real'
N_c=subs(n,s,sqrt(-1)*w);
P_ne=real(N_c);
P_no=imag(N_c);
D_c=subs(d,s,sqrt(-1)*w);
P_de=real(D_c);
P_do=imag(D_c);

g=P_ne*P_de+P_no*P_do
% 
%  syms epsilon
%  prog=sosprogram(w);
%  prog=sosdecvar(prog,[epsilon k_p k_i]);
%       
%       
%       prog=sosineq(prog,c1);
%       prog=sosineq(prog,c2);
%       prog=sosineq(prog,c3);
%       prog=sosineq(prog,c4);
% 
%       
%      prog=sosineq(prog,g-epsilon);
%      prog=sossetobj(prog,-epsilon);
%      prog=sossolve(prog);
%      
%      e_value=sosgetsol(prog, epsilon)
%      k_p1=sosgetsol(prog,k_p)
%      k_i1=sosgetsol(prog,k_i)

%  

%% OFP INDEX MAXIMIZED

syms xi;

    prog1=sosprogram(w);
    prog1=sosdecvar(prog1,[xi k_p k_i]);
    M=P_ne*P_de+P_no*P_do-xi*P_ne^2-xi*P_no^2;

    prog1=sosineq(prog1,c1);
    prog1=sosineq(prog1,c2);
    prog1=sosineq(prog1,c3);
    prog1=sosineq(prog1,c4);
     prog1=sosmatrixineq(prog1,M);
     prog1=sossetobj(prog1,-xi);
     prog1=sossolve(prog1);
     

     xi=double(sosgetsol(prog1,xi))
     k_p1=double(sosgetsol(prog1,k_p))
     k_i1=double(sosgetsol(prog1,k_i))

% optimal solution xi=0.5418 k_p1=1,k_p2=1

%% IFP INDEX MAXIMIZED
syms nu e;
nu_upper=1;
nu_lower=0;
error=1;


while(error>1e-4)
    nu=1/2*(nu_upper+nu_lower);
    prog1=sosprogram(w);
    prog1=sosdecvar(prog1,[e k_p k_i]);
    M=[1*P_ne*P_de+1*P_no*P_do sqrt(nu)*P_de sqrt(nu)*P_do;
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




    