function [Ainv] = inv_3by3(A)
% see mathworld.wolfram.com/MatrixInverse.html

C11=det([A(2,2),A(2,3);A(3,2),A(3,3)]);
C12=det([A(1,3),A(1,2);A(3,3),A(3,2)]);
C13=det([A(1,2),A(1,3);A(2,2),A(2,3)]);

C21=det([A(2,3),A(2,1);A(3,3),A(3,1)]);
C22=det([A(1,1),A(1,3);A(3,1),A(3,3)]);
C23=det([A(1,3),A(1,1);A(2,3),A(2,1)]);

C31=det([A(2,1),A(2,2);A(3,1),A(3,2)]);
C32=det([A(1,2),A(1,1);A(3,2),A(3,1)]);
C33=det([A(1,1),A(1,2);A(2,1),A(2,2)]);

Ainv=1/det(A)*[C11,C12,C13;C21,C22,C23;C31,C32,C33];
