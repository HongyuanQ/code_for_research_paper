function initial_displacement = find_initial_condition
global H M L1 L2 n1 n2 u
u = 1;
% Drillpipe specification
E1 = 210*10^9;      % E is modulus of elasticity, 210*10^9 N/m^2
A1 = 5.58*10^(-3);  % A is cross sectional area, 5.58*10^(-3) m^2 
L1 = 1000;            % L is the length of the pipe, 600 m
n1 = 21;            % n1 is the node of the pipe
l1 = L1/(n1-1);      % l is the length of the finite element of the pipe
I1 = 8.77*10^(-6);  % I is the moment of inertia, 8.77*10^(-6) m^4
G1 = 7.6923*10^10;  % G is the shear modulus, 7.6923*10^10 N/m^2
J1 = 1.754*10^(-5); % J is the polar moment of inertia for drillpipe, 1.754*10^(-5) m^4

Kpipe = zeros(2*(n1-1)+6,2*(n1-1)+6); % K is the global stiffness matrix

k1 = [E1*A1/l1       0         -E1*A1/l1        0;         % k1 is the local stiffness matrix
    
      0        G1*J1/l1         0         -G1*J1/l1;      
      
    -E1*A1/l1       0          E1*A1/l1        0;      
      
      0       -G1*J1/l1         0          G1*J1/l1] ;
  
k2 = [E1*A1/l1       0         -E1*A1/l1        0        0        0        0        0;         % k2 is the joint stiffness matrix
    
      0        G1*J1/l1         0            0        0     -G1*J1/l1      0        0;      
      
     -E1*A1/l1       0          E1*A1/l1        0        0        0        0        0;
     
      0           0          0            0        0        0        0        0;
      
      0           0          0            0        0        0        0        0
      
      0       -G1*J1/l1         0            0        0      G1*J1/l1      0        0
      
      0           0          0            0        0        0        0        0
      
      0           0          0            0        0        0        0        0];
  
  
 for i1=1:(n1-1)    % assemble stiffness matrix n-1 times
     if i1<=n1-2
         Kpipe((2*(i1-1)+1):(2*(i1+1)),(2*(i1-1)+1):(2*(i1+1))) = Kpipe((2*(i1-1)+1):(2*(i1+1)),(2*(i1-1)+1):(2*(i1+1)))+k1(:,:);
     end
     if i1 == n1-1
         Kpipe((2*(i1-1)+1):(2*(n1-1)+6),(2*(i1-1)+1):(2*(n1-1)+6)) = Kpipe((2*(i1-1)+1):(2*(n1-1)+6),(2*(i1-1)+1):(2*(n1-1)+6))+k2(:,:);
     end     
 end
 
 % mass matrix assemble
p = 7850;          % p is the density of the materia, 7850 kg/m^3

m1 =[p*A1*l1/3       0        p*A1*l1/6        0;         % m is the local mass matrix
    
      0         p*J1*l1/3        0        p*J1*l1/6;
      
     p*A1*l1/6       0        p*A1*l1/3        0;      
      
      0         p*J1*l1/6        0        p*J1*l1/3];
 
m2 =[p*A1*l1/3       0        p*A1*l1/6        0          0          0           0          0;         % m is the joint mass matrix
    
      0         p*J1*l1/3        0           0          0       p*J1*l1/6        0          0;
      
     p*A1*l1/6       0        p*A1*l1/3        0          0          0           0          0;
     
      0            0           0           0          0          0           0          0;
      
      0            0           0           0          0          0           0          0;
      
      0         p*J1*l1/6        0           0          0       p*J1*l1/3        0          0;
      
      0            0           0           0          0          0           0          0
      
      0            0           0           0          0          0           0          0];  
  
Mpipe = zeros(2*(n1-1)+6,2*(n1-1)+6); % M is the global mass matrix  
 for i1=1:(n1-1)    % assemble stiffness matrix n-1 times
     if i1<=n1-2
         Mpipe((2*(i1-1)+1):(2*(i1+1)),(2*(i1-1)+1):(2*(i1+1))) = Mpipe((2*(i1-1)+1):(2*(i1+1)),(2*(i1-1)+1):(2*(i1+1)))+m1(:,:);
     end
     if i1 == n1-1
         Mpipe((2*(i1-1)+1):(2*(n1-1)+6),(2*(i1-1)+1):(2*(n1-1)+6)) = Mpipe((2*(i1-1)+1):(2*(n1-1)+6),(2*(i1-1)+1):(2*(n1-1)+6))+m2(:,:);
     end     
 end

% Drillcollar specification
E = 210*10^9;      % E is modulus of elasticity, 210*10^9 N/m^2
A = 3.6483*10^(-2);% A is cross sectional area, 3.6483*10^(-2) 
L2 = 200;            % L is the length of the Drillcollar, 20 m
n2 = 11;            % n2 is the node of the system
l = L2/(n2-1);       % l is the length of the finite element of the Drillcollar
I = 1.324*10^(-4); % I is the moment of inertia, 1.324*10^(-4) m^4
G = 7.6923*10^10;  % G is the shear modulus, 7.6923*10^10 N/m^2
J = 2.648*10^(-4); % J is the polar moment of inertia for Drillcollar, 2.648*10^(-5) m^4

Kcollar = zeros(6*n2,6*n2); % K is the global stiffness matrix

k = [E*A/l         0           0           0           0          0       -E*A/l        0         0        0          0            0;         % k is the local stiffness matrix
      0       12*E*I/l^3       0           0           0      6*E*I/l^2     0       -12*E*I/l^3   0        0          0         6*E*I/l^2;
      0            0         12*E*I/l^3    0        -6*E*I/l^2    0         0           0     -12*E*I/l^3  0       -6*E*I/l^2      0;
      0            0           0         G*J/l         0          0         0           0         0      -G*J/l       0            0;
      0            0         -6*E*I/l^2    0         4*E*I/l      0         0           0       6*E*I/l^2  0        2*E*I/l        0;
      0        6*E*I/l^2       0           0           0      4*E*I/l       0        -6*E*I/l^2   0        0          0         2*E*I/l;
    -E*A/l         0           0           0           0          0        E*A/l        0         0        0          0            0;
      0      -12*E*I/l^3       0           0           0     -6*E*I/l^2     0        12*E*I/l^3   0        0          0        -6*E*I/l^2;
      0            0        -12*E*I/l^3    0         6*E*I/l^2    0         0           0      12*E*I/l^3  0        6*E*I/l^2      0;
      0            0           0        -G*J/l         0          0         0           0         0       G*J/l       0            0;
      0            0         -6*E*I/l^2    0         2*E*I/l      0         0           0       6*E*I/l^2  0        4*E*I/l        0;
      0        6*E*I/l^2       0           0           0      2*E*I/l       0        -6*E*I/l^2   0        0          0         4*E*I/l];
 
  
 for i1=1:(n2-1)    % assemble stiffness matrix n-1 times
     Kcollar((6*(i1-1)+1):(6*(i1+1)),(6*(i1-1)+1):(6*(i1+1))) = Kcollar((6*(i1-1)+1):(6*(i1+1)),(6*(i1-1)+1):(6*(i1+1)))+k(:,:);
 end

% mass matrix assemble
p = 7850;          % p is the density of the materia, 7850 kg/m^3
b = p*A*l/420;

m = [p*A*l/3       0           0           0           0          0       p*A*l/6       0         0        0          0            0;         % m is the local mass matrix
      0          156*b         0           0           0       22*l*b       0          54*b       0        0          0         -13*l*b;
      0            0         156*b         0        -22*l*b       0         0           0        54*b      0       13*l*b          0;
      0            0           0       p*J*l/3         0          0         0           0         0     p*J*l/6       0            0;
      0            0       -22*l*b         0         4*l^2*b      0         0           0       -13*l*b    0      -3*l^2*b         0;
      0         22*l*b         0           0           0      4*l^2*b       0         13*l*b      0        0          0         -3*l^2*b;
     p*A*l/6       0           0           0           0          0       p*A*l/3       0         0        0          0            0;
      0          54*b          0           0           0       13*l*b       0         156*b       0        0          0         -22*l*b;
      0            0         54*b          0        -13*l*b       0         0           0       156*b      0        22*l*b         0;
      0            0           0        p*J*l/6        0          0         0           0         0     p*J*l/3       0            0;
      0            0        13*l*b         0        -3*l^2*b      0         0           0       22*l*b     0        4*l^2*b        0;
      0         -13*l*b        0           0           0      -3*l^2*b      0        -22*l*b      0        0          0         4*l^2*b];
  
Mcollar = zeros(6*n2,6*n2); % M is the global mass matrix  
 for i1=1:(n2-1)    % assemble mass matrix n-1 times
    Mcollar((6*(i1-1)+1):(6*(i1+1)),(6*(i1-1)+1):(6*(i1+1))) = Mcollar((6*(i1-1)+1):(6*(i1+1)),(6*(i1-1)+1):(6*(i1+1)))+m(:,:);
 end 
 
K = zeros(2*(n1-1)+6*n2,2*(n1-1)+6*n2);
 for i1=1:2    % assemble global stiffness matrix
     if i1 == 1
       K((1:2*(n1-1)+6),(1:2*(n1-1)+6)) = Kpipe(:,:);
     end
     if i1 == 2
       K(((2*(n1-1)+1):(2*(n1-1)+6*n2)),((2*(n1-1)+1):(2*(n1-1)+6*n2))) = K(((2*(n1-1)+1):(2*(n1-1)+6*n2)),((2*(n1-1)+1):(2*(n1-1)+6*n2)))+Kcollar(:,:);
     end      
 end
 
 M = zeros(2*(n1-1)+6*n2,2*(n1-1)+6*n2);
 for i1=1:2    % assemble global mass matrix
     if i1 == 1
       M((1:2*(n1-1)+6),(1:2*(n1-1)+6)) = Mpipe(:,:);
     end
     if i1 == 2
       M(((2*(n1-1)+1):(2*(n1-1)+6*n2)),((2*(n1-1)+1):(2*(n1-1)+6*n2))) = M(((2*(n1-1)+1):(2*(n1-1)+6*n2)),((2*(n1-1)+1):(2*(n1-1)+6*n2)))+Mcollar(:,:);
     end      
 end 
 
C = 0.01*M+0.01*K;      % damping
 
 %************************************************************************************************************************************************% 
% consider the boundary condition
K((2*(n1-1)+2),(2*(n1-1)+2)) = K((2*(n1-1)+2),(2*(n1-1)+2))+10000000;   % stablizer 1
K((2*(n1-1)+3),(2*(n1-1)+3)) = K((2*(n1-1)+3),(2*(n1-1)+3))+10000000;
K((2*(n1-1)+5),(2*(n1-1)+5)) = K((2*(n1-1)+5),(2*(n1-1)+5))+10000000;
K((2*(n1-1)+6),(2*(n1-1)+6)) = K((2*(n1-1)+6),(2*(n1-1)+6))+10000000;

K((2*(n1-1)+6*(n2-1)+2),(2*(n1-1)+6*(n2-1)+2)) = K((2*(n1-1)+6*(n2-1)+2),(2*(n1-1)+6*(n2-1)+2))+10000000;       % stablizer 2
K((2*(n1-1)+6*(n2-1)+3),(2*(n1-1)+6*(n2-1)+3)) = K((2*(n1-1)+6*(n2-1)+3),(2*(n1-1)+6*(n2-1)+3))+10000000;
K((2*(n1-1)+6*(n2-1)+5),(2*(n1-1)+6*(n2-1)+5)) = K((2*(n1-1)+6*(n2-1)+5),(2*(n1-1)+6*(n2-1)+5))+10000000;
K((2*(n1-1)+6*(n2-1)+6),(2*(n1-1)+6*(n2-1)+6)) = K((2*(n1-1)+6*(n2-1)+6),(2*(n1-1)+6*(n2-1)+6))+10000000;
%************************************************************************************************************************************************%Kpipe

N = length(M);
 
%***********************************************************************%
% rock stiffness
E_Rock = 40*10^9;
nu = 0.25;
G_Rock = E_Rock/(2*(1+nu));
Kc = G_Rock*0.2/(1-nu);
K(N-5,N-5) = K(N-5,N-5)+Kc;

H = [zeros(N)     eye(N);
      -M\K         -M\C];
  
tspan = [0 40];
x0 = zeros(2*N,1);
options = [];

tic
[t,x] = ode23t('displacement_under_gravity',tspan,x0,options); 
toc

initial_displacement = x(length(t),1:N);

