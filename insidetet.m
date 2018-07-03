function [ckin] = insidetet(v1,v2,v3,v4,p)
% v1-v4 are the vertices of the tetrahedron
% p is the sample point
% ckin is a boolean variable: 1=inside, 0=outside

ovec=ones(4,1);

det0=det( [[v1';v2';v3';v4'],ovec] );
det1=det( [[p';v2';v3';v4'],ovec] );
det2=det( [[v1';p';v3';v4'],ovec] );
det3=det( [[v1';v2';p';v4'],ovec] );
det4=det( [[v1';v2';v3';p'],ovec] );

% sign(det0)
% sign(det1)
% sign(det2)
% sign(det3)
% sign(det4)
ssum=sign(det0)+sign(det1)+sign(det2)+sign(det3)+sign(det4);
ckin=abs(ssum)==5;