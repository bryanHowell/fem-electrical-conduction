function [delN] = dirshape_ele(pele)
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

% pre-allocation
% - 10 rows: one row per function
% - 4 columns: N, dNde, dNdn, dNdz
delN=zeros(3,10);

% shape function derivatives (w/r to psi, eta, and zeta)
delN(1,:)=[4*(e+n+c)-3,4*e-1,0,0,4*(1-2*e-n-c),-4*n,-4*c,4*n,0,4*c];
delN(2,:)=[4*(e+n+c)-3,0,4*n-1,0,-4*e,4*(1-e-2*n-c),-4*c,4*e,4*c,0];
delN(3,:)=[4*(e+n+c)-3,0,0,4*c-1,-4*e,-4*n,4*(1-e-n-2*c),0,4*n,4*e];
