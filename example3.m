clear all;
syms z;
D=z-2;
double(solve(D==0))
N=z;

syms k_p k_i k_d
C=k_p+k_i/(z-0.5);%+k_d*s;

c1=(1-k_p)*(k_p-0.1);
c2=(2-k_i)*(k_i-1);
%c3=(1-k_d)*(k_d+1);

G=(N/D)/(1+(N/D)*C);


[n,d] = numden(simplify(G));
[coeff,ps]=coeffs(d,z);


[J,JT]=modified_jury(coeff);
[a,b]=size(JT);

length=a;
fcolumn=sym('fcolumn',[length,1]);


for i=1:length     %the 1,1-th entry must be 1, ignore it
    fcolumn(i)=JT(i,1)
end

syms e
z=[k_p, k_i];

% check the posivity of fcolumn(1)&(2) manually
     prog=sosprogram(z);
     prog=sosdecvar(prog,e)
    [prog,s1]=sossosvar(prog,monomials(z,[0:1])); % choose the degree of lowerbound polynomial Xi
    [prog,s2]=sossosvar(prog,monomials(z,[0:1]));
    %[prog,s3]=sossosvar(prog,monomials(z,[0:0]));
     g=fcolumn(3)-c1*s1-c2*s2-e;
     prog=sosineq(prog,g);
     prog=sossetobj(prog,-e);
     prog=sossolve(prog)
     e_value=sosgetsol(prog,e)
     s1_value=simplify(sosgetsol(prog,s1))
     s2_value=sosgetsol(prog,s2)

 
  