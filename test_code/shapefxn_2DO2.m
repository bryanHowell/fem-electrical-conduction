function [Ndata] = shapefxn_2DO2(e,n)
% triangular elements in 2D
% - 6 nodes per element
% => six function evaluations per element

% element domain:
% e, psi [0,1]
% n, eta [0,1]

% pre-allocation
% - 6 rows: one row per function
% - 3 columns: N, dNde, dNdn
Ndata=zeros(6,3);

% shape function 
Ndata(:,1)=[(1-e-n)*(1-2*e-2*n);e*(2*e-1);n*(2*n-1);...
    4*e*(1-e-n);4*e*n;4*n*(1-e-n)];

% shape function derivatives (w/r to psi and eta)
Ndata(:,2)=[4*e+4*n-3;4*e-1;0;4*(1-2*e-n);4*n;-4*n];
Ndata(:,3)=[4*e+4*n-3;0;4*n-1;-4*e;4*e;4*(1-e-2*n)];



