function [Ndata] = shapefxn_3DO2(e,n,c)
% triangular elements in 3D
% - 10 nodes per element
% => ten function evaluations per element

% element domain:
% e, psi [1,0,0]
% n, eta [0,1,0]
% c, zeta [0,0,1]

% pre-allocation
% - 10 rows: one row per function
% - 4 columns: N, dNde, dNdn, dNdz
Ndata=zeros(10,4);

% shape function 
Ndata(:,1)=[(2*(1-e-n-c)-1)*(1-e-n-c);...
            (2*e-1)*e;(2*n-1)*n;(2*c-1)*c;...
            4*(1-e-n-c)*e;4*(1-e-n-c)*n;4*(1-e-n-c)*c;...
            4*e*n;4*n*c;4*e*c];

% shape function derivatives (w/r to psi, eta, and zeta)
Ndata(:,2)=[4*(e+n+c)-3;4*e-1;0;0;4*(1-2*e-n-c);-4*n;-4*c;4*n;0;4*c];
Ndata(:,3)=[4*(e+n+c)-3;0;4*n-1;0;-4*e;4*(1-e-2*n-c);-4*c;4*e;4*c;0];
Ndata(:,4)=[4*(e+n+c)-3;0;0;4*c-1;-4*e;-4*n;4*(1-e-n-2*c);0;4*n;4*e];
