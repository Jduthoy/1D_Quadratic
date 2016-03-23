function [ F ] = FEMquad( H,g,ne,ln,f )
%%%%%%% Julius Duthoy%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Math 520     %%%%%%%%%%%%%%%%%%%%%%%%%%
%Test 1: FEMquad(0,0,4,3,@(x) cos(x))
% u" + u = cos(x)
% B.C. -u'(0) = H=0, u(1) = g=0

%Test 2: FEMquad(2,3,4,3,@(x) cos(x))
%%%%%%%%%%List of Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%ln  =# Local Nodes
%gn  =# Global Nodes 
%ne  =# elements- input into the function
%n_eq=# Nonzero Global Basis Fncts

gn=2*ne+1;

%%%%%%%%%% Set up Global Nodes & Constant Interval h %%%%%%%%%%%%%%%%%%
x=[0:1/(2*ne):1];
h=x(3:2:gn)-x(1:2:gn-2);

%%%%%%%%%% ID(A)=P Matrix%%%%%%%%%%%%%%%%%%%%%%%%%% 
u(1)=g;
n_eq=0;%Determined from ID Matrix below

for A=1:gn
    if A==gn
        ID(A)=0; %only works for u(1)=g not u(0)=g
    else
        n_eq=n_eq+1;
        ID(A)=n_eq;
    end
end
ID;

%%%%%% IEN(a,e)=A Matrix %%%%%%%%%%%%%%%%%%%%%%%%
for e=1:ne 
    for a=1:ln
        IEN(a,e)=a+2*(e-1);
    end
end
IEN;

%%%%%%%%% LM(a,e)=ID(IEN) Matrix %%%%%%%%%%%%%%%%%%%%%%%%
LM=ID(IEN);


%%%%%%% Initialize Local Basis fncts %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Derived from Lagrange %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N1=@(z) 2*z.^2-3*z+1;
dN1=@(z) 4*z-3;

N2=@(z) -4*z.^2 + 4*z;
dN2=@(z) -8*z+4;

N3=@(z) 2*z.^2-z;
dN3=@(z) 4*z-1;


N={'2*z.^2-3*z+1','-4*z.^2 + 4*z','2*z.^2-z'};
dN={'4*z-3','-8*z+4','4*z-1'};

%%%%%%%%%Stiffness Matrix and F Vector%%%%%%%%%%%%%%%%%%%%%%%%%
K=zeros(n_eq,n_eq); 
F=zeros(n_eq,1);

for e=1:ne
       %{
       k(1,1)=-(1/h(e)).*quad(@(z) dN1(z).*dN1(z), 0, 1)+h(e).*quad(@(z) N1(z).*N1(z), 0, 1);
       k(1,2)=-(1/h(e)).*quad(@(z) dN1(z).*dN2(z), 0, 1)+h(e).*quad(@(z) N1(z).*N2(z), 0, 1);
       k(1,3)=-(1/h(e)).*quad(@(z) dN1(z).*dN3(z), 0, 1)+h(e).*quad(@(z) N1(z).*N3(z), 0, 1);
       k(2,1)=-(1/h(e)).*quad(@(z) dN2(z).*dN1(z), 0, 1)+h(e).*quad(@(z) N2(z).*N1(z), 0, 1);
       k(2,2)=-(1/h(e)).*quad(@(z) dN2(z).*dN2(z), 0, 1)+h(e).*quad(@(z) N2(z).*N2(z), 0, 1);
       k(2,3)=-(1/h(e)).*quad(@(z) dN2(z).*dN3(z), 0, 1)+h(e).*quad(@(z) N2(z).*N3(z), 0, 1);
       k(3,1)=-(1/h(e)).*quad(@(z) dN3(z).*dN1(z), 0, 1)+h(e).*quad(@(z) N3(z).*N1(z), 0, 1);
       k(3,2)=-(1/h(e)).*quad(@(z) dN3(z).*dN2(z), 0, 1)+h(e).*quad(@(z) N3(z).*N2(z), 0, 1);
       k(3,3)=-(1/h(e)).*quad(@(z) dN3(z).*dN3(z), 0, 1)+h(e).*quad(@(z) N3(z).*N3(z), 0, 1);
    %}
       %Set Up Big K
       for a=1:ln
            if LM(a,e)~=0
             fnct=inline(N{1,a});
             F(LM(a,e))=F(LM(a,e))+h(e).*quad(@(z) f(x(LM(1,e))+h(e).*z).*fnct(z), 0, 1);
                for b=1:ln
                    if LM(b,e)~=0
                        fnct_a=inline(N{1,a});
                        fnct_b=inline(N{1,b});
                        der_a=inline(dN{1,a});
                        der_b=inline(dN{1,b});
                        k(a,b)=-(1/h(e)).*quad(@(z) der_a(z).*der_b(z), 0, 1)+h(e).*quad(@(z) fnct_a(z).*fnct_b(z), 0, 1);
                        K(LM(a,e),LM(b,e))=K(LM(a,e),LM(b,e))+k(a,b);
                    end
                end
            end
       end    

       % F matrix
%{
       %if LM(1,e)~=0
            F(LM(1,e))=F(LM(1,e))+h(e).*quad(@(z) f(x(LM(1,e))+h(e).*z).*(2*z.^2-3*z+1), 0, 1);
       %end

       %if LM(2,e)~=0
            F(LM(2,e))=F(LM(2,e))+h(e).*quad(@(z) f(x(LM(1,e))+h(e).*z).*(-4*z.^2+4*z), 0, 1);
       %end

       if LM(3,e)~=0
            F(LM(3,e))=F(LM(3,e))+h(e).*quad(@(z) f(x(LM(1,e))+h(e).*z).*(2*z.^2-z), 0, 1);
       end
%}
      % Placing B.C.s in proper entries
       if e==ne
            F(LM(1,e))=F(LM(1,e))-k(1,3)*g; 
            F(LM(2,e))=F(LM(2,e))-k(2,3)*g;
            %F(LM(3,e))=F(LM(3,e))-k(3,3)*g;
       end

       if e==1
            F(LM(1,e))=F(LM(1,e))-H*N1(0);
       end
end

% Solve for d
d = K\F;

% Get u from d
for A=1:gn
    if ID(A)~=0
        u(A)=d(ID(A));
    else
        u(A)=g;
    end
end

%%%%%%%%%%Comparison to Exact Soln%%%%%%%%%%%%%%%%%%%%%%%%%
%sol=@(z) (1/2)*(z*sin(z)-tan(1)*cos(z));
sol=@(z) (1/2)*[(z-4)*sin(z)+3*(2+sin(1))*sec(1)*cos(z)];
UU=arrayfun(sol, x);%UU=Actual u calculated by hand
%arrayfun(function,inputs) uses x array as inputs to sol to form the actual

%%%%%%%%%% Plots the Actual U and my FEM approximation
plot(x, u, 'm*', x, UU, 'k')
legend('Approx u', 'Actual u')
str=sprintf('u"+u=cos(x), u(1)=%f, -u_x(0)=%f',g,H);
% Couldn't get f to change to whatever the user inputs
title(str);
xlabel('x')
ylabel('u')

%format long
e=max(abs(UU-u)) %Calculates the error

%K
end

