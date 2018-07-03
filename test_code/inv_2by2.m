function [Ainv] = inv_2by2(A)

a=A(1,1);
b=A(1,2);
c=A(2,1);
d=A(2,2);
Ainv=1/(a*d-b*c)*[d,-b;-c,a];