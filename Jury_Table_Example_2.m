clear all;
syms x1 x2 k s;
x=[x1 x2];
syms e; 

Ap=[1.1*x1+1.2*x2];
Bp=[5*x1+4*x2];
Cp=[-1-0.5*x2];
 

Dc=[k];


E=[Ap];
F=[Bp];
G=[Dc*Cp];


PHI=kron((E+F*II*G),(E+F*II*G))+S*kron((F*G),(F*G));
chapoly=det(s*eye(1)-PHI);
[coeff,ps]=coeffs(chapoly, s);
%coeff=fliplr(coeff);
[J,JuryTable]=jury(coeff);
[a,b]=size(JuryTable);
length=(a+1)/2;
syms X;
fcolumn=X*ones(length-1,1);
for i=1:(length-1)     %the 1,1-th entry must be 1, ignore it
    fcolumn(i)=JuryTable(2*i+1,1);
end
for i=1:length-1

  fcolumn(i)=numden(fcolumn(i));
  %fcolumn(i)=homogenneous(fcolumn(i),x);% the column of numerators need to be checked
end

z=[x.'; k].';
gamma=sym('gamma',[length-1,1]);
g=sym('g',[length-1,1]);
gm=sym('gm',[length-1,1]);
gammam=sym('gammam',[length-1,1]);

Fi=1;
c1=1-k^2;

m=2;
while(m<=8) %degree of the scalar polynomial variable gamma


    prog=sosprogram(z);
    [prog,Xi]=sospolyvar(prog,monomials(k,0:4),'wscoeff');
    %prog=sosdecvar(prog,[A;B;C;D]);
   
 
   
    H=@(k)Fi*Xi;
     for i=1:(length-1)
       [prog,gamma(i)]=sospolyvar(prog,monomials([x1 x2 k^2],[0:m]),'wscoeff');  
       %   [prog, Xi]=sospolyvar(prog,monomials(k,[0:3]),'wscoeff');
        g(i)=fcolumn(i)-Xi-c1*gamma(i);
        
        gm(i)=homogeneous(g(i),x);
           
        if m==0
           gammam(i)=gamma(i);
            else
             gammam(i)=homogeneous(gamma(i),x); 
        end
        
        gm(i)=subs(gm(i),[x1 x2],[x1^2 x2^2]);
        gammam(i)=subs(gammam(i),[x1 x2],[x1^2 x2^2]);
     
        prog=sosineq(prog,gammam(i),'sparse');    
        prog=sosineq(prog,gm(i),'sparse');
   
        
     end
    
    prog=sossetobj(prog,-int(H,k,-1,1));
    prog=sossolve(prog)
    
       if (prog.solinfo.info.pinf==0 && prog.solinfo.info.numerr~=2)
       Xiv=sosgetsol(prog,Xi)
       %if subs(Xiv,k,0)>0 && subs(Xiv,k,1)>0 && subs(Xiv,k,-1)>0
          
           Xi=Xiv;
           handle=matlabFunction(-Xi);
           [K,v]=fminbnd(handle,-1,1)
             if v<0
             disp('found the positive lower bound with the i-th entry')  
             break;
             end
       
      end
  
  m=m+1;
 end