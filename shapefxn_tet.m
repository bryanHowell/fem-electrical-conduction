function [N] = shapefxn_tet(pele)
% triangular elements in 3D
% - 10 nodes per element
% => ten function evaluations per element

% pele, point in parent/element domain
% e, psi
% n, eta
% c, zeta

e=pele(1);
n=pele(2);
c=pele(3);

N=[(2*(1-e-n-c)-1)*(1-e-n-c),...
   (2*e-1)*e,(2*n-1)*n,(2*c-1)*c,...
    4*(1-e-n-c)*e,4*(1-e-n-c)*n,4*(1-e-n-c)*c,...
    4*e*n,4*n*c,4*e*c];
