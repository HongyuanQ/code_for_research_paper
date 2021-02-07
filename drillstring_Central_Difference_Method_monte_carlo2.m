% this code will use central difference method to simulate dynamic response of the drillstring 

clear all;
close all;
clc

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
      
      0       -G1*J1/l1        0          G1*J1/l1] ;
  
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
L2 = 200;            % L is the length of the Drillcollar, 200 m
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

M_stabilizer = 500;
M((2*(n1-1)+1),(2*(n1-1)+1)) = M((2*(n1-1)+1),(2*(n1-1)+1))+M_stabilizer;
M((2*(n1-1)+6*(n2-1)+1),(2*(n1-1)+6*(n2-1)+1)) = M((2*(n1-1)+6*(n2-1)+1),(2*(n1-1)+6*(n2-1)+1))+M_stabilizer;
r_inb = 0.03;

C = 0.01*M+0.01*K;      % damping

%************************************************************************************************************************************************%Kpipe

N = length(M);
delt = 2*pi/10/sqrt(max(eig(K,M)))

initial_displacement = find_initial_condition;
displacement_X = zeros(N,2);
displacement_X(:,2) = initial_displacement';
velosity_X = zeros(N,1); 
intial1 = displacement_X;

wd = 15; % this is the rotary table speed, stick-slip 15, no stick-slip 30

for i1 = 1:(n1+n2-1)
    if i1 <= n1-1
        velosity_X(2*(i1-1)+2) = wd;
    end
    if i1 > n1-1
        velosity_X(2*(n1-1)+6*(i1-n1)+4) = wd;
    end    
end
intial2 = velosity_X;

% assemble constant force without DC motor
A1 = 5.58*10^(-3);  % A1 is pipe sectional area, 5.58*10^(-3) m^2
A2 = 3.6483*10^(-2);% A2 is collar sectional area, 3.6483*10^(-2)
p = 7850;           % p is the density of the materia, 7850 kg/m^3

l1 = L1/(n1-1);
l2 = L2/(n2-1);
f = zeros(2*(n1-1)+6*n2,1);

% constant force
for i1 = 1:(n1+n2-2)
    if i1 <=(n1-1)
        f([(2*(i1-1)+1) (2*i1+1)]) = f([(2*(i1-1)+1) (2*i1+1)])+[p*l1*A1*9.84/2;p*l1*A1*9.84/2];
    end
    if i1 >(n1-1)
        f([(2*(n1-1)+6*(i1-n1)+1) (2*(n1-1)+6*(i1-n1+1)+1)]) = f([(2*(n1-1)+6*(i1-n1)+1) (2*(n1-1)+6*(i1-n1+1)+1)])+[p*l2*A2*9.84/2;p*l2*A2*9.84/2];
    end    
end

f2 = f;
Fh   = 0.6*(p*L1*A1*9.84+p*L2*A2*9.84+M_stabilizer*2*9.84);  % hookload
f(1) = f(1)-Fh;
f(2*(n1-1)+1) = f(2*(n1-1)+1)+M_stabilizer*9.84;
f(2*(n1-1)+6*(n2-1)+1) = f(2*(n1-1)+6*(n2-1)+1)+M_stabilizer*9.84;
f2(1) = f2(1) - 0.8*(p*L1*A1*9.84+p*L2*A2*9.84+M_stabilizer*2*9.84);
f2(2*(n1-1)+1) = f2(2*(n1-1)+1)+M_stabilizer*9.84;
f2(2*(n1-1)+6*(n2-1)+1) = f2(2*(n1-1)+6*(n2-1)+1)+M_stabilizer*9.84;

% monte carlo simulation and central difference method
dt = 0.0002;
N1 = M+1/2*dt*C;
N2 = N1\(2*M-dt^2*K);
N3 = N1\(1/2*dt*C-M);

S0 = 0;       % S0 should be larger than 0 for random excitation 
Sample = 1;   % Sample should be larger than 1 for random excitation
t = 0:dt:80;
N4 = length(t);
SampleX = zeros(Sample,N4,6);
SampleX(:,1,1) = velosity_X(N-2)*ones(Sample,1);                     % bit rotation
SampleX(:,1,2) = displacement_X(N-5,2)*ones(Sample,1);               % bit axial displacement
SampleX(:,1,3) = displacement_X(N-2,2)*ones(Sample,1);               % bit rotational displacement
SampleX(:,1,4) = velosity_X(N-5)*ones(Sample,1);                     % bit axial velosity
SampleX(:,1,5) = displacement_X(N-4,2)*ones(Sample,1);               % bit lateral_y displacement
SampleX(:,1,6) = displacement_X(N-3,2)*ones(Sample,1);               % bit lateral_z displacement
%***********************************************************************%
% rock stiffness
E_Rock = 40*10^9;
nu = 0.25;
G_Rock = E_Rock/(2*(1+nu));
Kc = G_Rock*0.2/(1-nu);

% rock damping
rho = 2100;
c = 3.4*0.2^2*sqrt(G_Rock*rho)/(1-nu);

f0 = p*L2*A2*9.84+p*L1*A1*9.84-Fh+M_stabilizer*2*9.84;
average_speed = wd;
ROP = 1.35*(10^(-8))*f0*sqrt(average_speed)-1.9*10^-4;
zetac = 2*pi*ROP/average_speed;
rb = 0.22;
sigma = 1;
ConPa = sigma*sqrt(zetac/rb);          % constant parameter in TOB equation
%***********************************************************************%

%***********************************************************************%
%DC motor constant
Km = 6;
n = 7.2;
L = 0.005;
Rm = 0.01;
desired_speed = wd;
Vc = desired_speed*Km*n;
SampleI = zeros(1,3);
SampleI(1) = 500;
SampleI(2) = 500;
T_0 = Km*n*SampleI(2);
initial3 = SampleI;

B1 = zeros(N,1);
B2 = ones(n2,1);
B1(2*(n1-1)+2:6:N) = B2;

B9 = zeros(N,1);
B10 = ones(n2-1,1);
B9(2*(n1-1)+2:6:N-6) = B10;
B11 = zeros(n2-1,N);
for i1 = 1:n2-1
    B11(i1,2*(n1-1)+6*(i1-1)+2) = 1;
end

B12 = zeros(N,n2-1);
for i1 = 1:n2-1
    B12(2*(n1-1)+6*(i1-1)+2,i1) = 1;
end

B14 = zeros(N,n2-1);
for i2 = 1:n2-1
    B14(2*(n1-1)+6*(i2-1)+4,i2) = 1;
end

B15 = zeros(N,1);
B16 = ones(n2-1,1);
B15(2*(n1-1)+4:6:N-6) = B16;

force = 0; % this force is for lateral excitation on the bit

N5 = floor(N4/2);
tic
for i1 = 1:Sample
    displacement_X = intial1;
    velosity_X = intial2;
    SampleI = initial3;
        
    % 0.3 is the friction coefficient, 0.001 is s0 in eq23 of the paper 
    F0 = f+[zeros((N-3),1);1;0;0]*randn()*sqrt(2*pi*S0/dt)-(rb*(Kc*(displacement_X(N-5,2)-0.001*sin(displacement_X(N-2,2)))+c*(velosity_X(N-5)-0.001*velosity_X(N-2)*cos(displacement_X(N-2,2))))*...
         (displacement_X(N-5,2)-0.001*sin(displacement_X(N-2,2))>=0)*(0.06*(tanh(velosity_X(N-2))+(2*velosity_X(N-2))/(1+(velosity_X(N-2))^2)+0.01*velosity_X(N-2))+ConPa))...
         *[zeros((N-3),1);1;0;0]-(Kc*(displacement_X(N-5,2)-0.001*sin(displacement_X(N-2,2)))+c*(velosity_X(N-5)-0.001*velosity_X(N-2)*cos(displacement_X(N-2,2))))*...
         (displacement_X(N-5,2)-0.001*sin(displacement_X(N-2,2))>=0)*[zeros((N-6),1);1;0;0;0;0;0]+[0;1;zeros(N-2,1)]*T_0+...
         [zeros((N-4),1);1;0;0;0]*M_stabilizer*velosity_X(N-2)^2*r_inb*sin(displacement_X(N-2,2))-[zeros((N-5),1);1;0;0;0;0]*M_stabilizer*velosity_X(N-2)^2*r_inb*cos(displacement_X(N-2,2))-...
         [zeros((N-3),1);1;0;0]*(M_stabilizer*velosity_X(N-2)^2*r_inb+force*sin(t(1)))*0.3*0.1078*sign(velosity_X(N-2))-...
         Kc*(B1.*displacement_X(:,2)-0.05).*(B1.*displacement_X(:,2)>0.05)-Kc*(B1.*displacement_X(:,2)+0.05).*(B1.*displacement_X(:,2)<-0.05)-...
         0.3*0.1078*B14*((Kc*(B11*(B9.*displacement_X(:,2))-0.05)).*(B11*(B9.*displacement_X(:,2))>0.05)).*sign(B15.*velosity_X)-...
         0.3*0.1078*B14*((Kc*(B11*(B9.*displacement_X(:,2))+0.05)).*(B11*(B9.*displacement_X(:,2))<-0.05)).*sign(B15.*velosity_X)+...
         [zeros((N-5),1);1;0;0;0;0]*force*sin(t(1));         
     
    acc_X = M\(F0-C*velosity_X-K*displacement_X(:,2));
    displacement_X(:,1) = displacement_X(:,2)-dt*velosity_X+dt^2/2*acc_X;
    for i2 = 2:N4
        next_X = dt^2*(N1\F0)+N2*displacement_X(:,2)+N3*displacement_X(:,1);
        
        rotary_table_velosity = velosity_X(2);
        
        velosity_X = velosity_X+acc_X*dt;
        displacement_X(:,1) = displacement_X(:,2);
        displacement_X(:,2) = next_X;
        
        SampleX(i1,i2,1) = velosity_X(N-2);
        SampleX(i1,i2,2) = displacement_X(N-5,2);
        SampleX(i1,i2,3) = displacement_X(N-2,2);
        SampleX(i1,i2,4) = velosity_X(N-5);
        SampleX(i1,i2,5) = displacement_X(N-4,2);
        SampleX(i1,i2,6) = displacement_X(N-3,2);
        
        % -- increase desired table speed to eliminate stick-slip in the middle of simulation --%
        %if i2>N5
            %Vc = 30*Km*n;
        %end
        % ********************************************************************************
        
        % -- increase hookload to eliminate stick-slip in the middle of simulation --% 
        %if i2>N5
            %f = f2;
        %end
        % ****************************************************************************   

        SampleI(3) = 1/L*(2*dt*Vc-2*dt*Km*n*rotary_table_velosity-2*Rm*SampleI(2)*dt+L*SampleI(1));
        SampleI(1) = SampleI(2);
        SampleI(2) = SampleI(3);
        T = Km*n*SampleI(2);
        
                
        F0 = f+[zeros((N-3),1);1;0;0]*randn()*sqrt(2*pi*S0/dt)-(rb*(Kc*(displacement_X(N-5,2)-0.001*sin(displacement_X(N-2,2)))+c*(velosity_X(N-5)-0.001*velosity_X(N-2)*cos(displacement_X(N-2,2))))*...
         (displacement_X(N-5,2)-0.001*sin(displacement_X(N-2,2))>=0)*(0.06*(tanh(velosity_X(N-2))+(2*velosity_X(N-2))/(1+(velosity_X(N-2))^2)+0.01*velosity_X(N-2))+ConPa))...
         *[zeros((N-3),1);1;0;0]-(Kc*(displacement_X(N-5,2)-0.001*sin(displacement_X(N-2,2)))+c*(velosity_X(N-5)-0.001*velosity_X(N-2)*cos(displacement_X(N-2,2))))*...
         (displacement_X(N-5,2)-0.001*sin(displacement_X(N-2,2))>=0)*[zeros((N-6),1);1;0;0;0;0;0]+[0;1;zeros(N-2,1)]*T+...
         [zeros((N-4),1);1;0;0;0]*M_stabilizer*velosity_X(N-2)^2*r_inb*sin(displacement_X(N-2,2))-[zeros((N-5),1);1;0;0;0;0]*M_stabilizer*velosity_X(N-2)^2*r_inb*cos(displacement_X(N-2,2))-...
         [zeros((N-3),1);1;0;0]*(M_stabilizer*velosity_X(N-2)^2*r_inb+force*sin(t(i2)))*0.3*0.1078*sign(velosity_X(N-2))-...
         Kc*(B1.*displacement_X(:,2)-0.05).*(B1.*displacement_X(:,2)>0.05)-Kc*(B1.*displacement_X(:,2)+0.05).*(B1.*displacement_X(:,2)<-0.05)-...
         0.3*0.1078*B14*((Kc*(B11*(B9.*displacement_X(:,2))-0.05)).*(B11*(B9.*displacement_X(:,2))>0.05)).*sign(B15.*velosity_X)-...
         0.3*0.1078*B14*((Kc*(B11*(B9.*displacement_X(:,2))+0.05)).*(B11*(B9.*displacement_X(:,2))<-0.05)).*sign(B15.*velosity_X)+...
         [zeros((N-5),1);1;0;0;0;0]*force*sin(t(i2));
     
        acc_X = M\(F0-C*velosity_X-K*displacement_X(:,2));       
    end    
end
toc

tic
mean_X = mean(SampleX(:,:,2),1);
mean_V = mean(SampleX(:,:,1),1);
TOB_0 = (Kc*(SampleX(:,:,2)-0.001*sin(SampleX(:,:,3)))+c*(SampleX(:,:,4)-0.001*SampleX(:,:,1).*cos(SampleX(:,:,3)))).*...
     (SampleX(:,:,2)>=0.001*sin(SampleX(:,:,3))).*0.22.*(0.06*(tanh(SampleX(:,:,1))+2*SampleX(:,:,1)./(1+1*SampleX(:,:,1).^2)...
      +0.01*SampleX(:,:,1))+1*sqrt(2*pi*(1.35*10^-8*f0*sqrt(average_speed)-1.9*10^-4)./(average_speed*0.22)))+...
      (M_stabilizer*SampleX(:,:,1).^2*r_inb+repmat(force*sin(t),[Sample,1]))*0.3*0.1078.*sign(SampleX(:,:,1));
TOB = mean(TOB_0,1);
WOB_0 = (Kc*(SampleX(:,:,2)-0.001*sin(SampleX(:,:,3)))+c*(SampleX(:,:,4)-0.001*SampleX(:,:,1).*cos(SampleX(:,:,3)))).*(SampleX(:,:,2)>=0.001*sin(SampleX(:,:,3)));
WOB = mean(WOB_0,1);
toc

figure(1)
subplot(3,1,1)
plot(t,mean_V,'r-')
xlabel('Time,sec')
ylabel('Speed (rad/s)')
legend('Bit Speed')
grid
subplot(3,1,2)
plot(t,TOB,'r-')
xlabel('Time,sec')
ylabel('Torque (Nm)')
legend('Bit torque')
grid
subplot(3,1,3)
plot(t,WOB,'r-')
xlabel('Time,sec')
ylabel('Weight on Bit (N)')
grid 

figure(2)
plot(t,mean_V,'k-')
xlabel('Time,sec')
ylabel('Mean Bit Speed (rad/s)')
legend('Mean Bit Speed')
grid
axes('position',[0.55,0.2,0.3,0.25]);
plot(t(70/0.0002+1:N4),mean_V(70/0.0002+1:N4),'k-')



% ----- only used for random excitations
% variance_w = std(SampleX(:,:,1),0,1);
% variance_x = std(SampleX(:,:,2),0,1);
% variance_y = std(SampleX(:,:,5),0,1);
% variance_z = std(SampleX(:,:,6),0,1);
% 
% figure(3)
% plot(t,variance_w,'k-')
% xlabel('Time,sec')
% ylabel('Bit Speed Standard deviation (rad/s)')
% legend('Bit Speed Standard deviation')
% grid
% axes('position',[0.55,0.2,0.3,0.25]);
% plot(t(70/0.0002+1:N4),variance_w(70/0.0002+1:N4),'k-')
% 
% figure(4)
% plot(t,variance_x,'k-')
% xlabel('Time,sec')
% ylabel('Bit displacement axial Standard deviation (m)')
% legend('Bit displacement axial Standard deviation')
% grid
% axes('position',[0.55,0.55,0.3,0.25]);
% plot(t(70/0.0002+1:N4),variance_x(70/0.0002+1:N4),'k-')
% 
% figure(5)
% plot(t,variance_y,'k-')
% xlabel('Time,sec')
% ylabel('Bit displacement lateral Standard deviation (m)')
% legend('Bit displacement lateral Standard deviation (y)')
% grid
% axes('position',[0.55,0.55,0.3,0.25]);
% plot(t(70/0.0002+1:N4),variance_y(70/0.0002+1:N4),'k-')
% 
% figure(6)
% plot(t,variance_z,'k-')
% xlabel('Time,sec')
% ylabel('Bit displacement lateral Standard deviation (m)')
% legend('Bit displacement lateral Standard deviation (z)')
% grid
% axes('position',[0.55,0.55,0.3,0.25]);
% plot(t(70/0.0002+1:N4),variance_z(70/0.0002+1:N4),'k-')
% ********************************************************************

figure(7)
plot(mean(SampleX(:,:,5),1),mean(SampleX(:,:,6),1),'k')
axis equal
xlabel('Y direction (m)')
ylabel('Z direction (m)')
legend('Bit trajectory')
grid



