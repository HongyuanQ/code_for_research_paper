function xprime = displacement_under_gravity(t,x)
global H M L1 L2 n1 n2 u
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

N = length(M);
Fh   = 0.6*(p*L1*A1*9.84+p*L2*A2*9.84);  % hookload
f(1) = f(1)-Fh;
Fc   = [zeros(N,1);M\f];

xprime = H*x+Fc*u;