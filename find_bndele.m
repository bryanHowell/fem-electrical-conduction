function [elein] = find_bndele(pcloud,bndele,rc,abc)

% pcloud, collection of points (# of points by 3)
% bndele, points (columns) comprising each element (row)
% rc, center of ellipsoid in 3D
% abc, xyz radii of ellipsoid, respectively

% see which points are within the specified ellipsoid
npts=size(pcloud,1);
ovec=ones(npts,1);
indx_in=( (pcloud(:,1)-rc(1)*ovec).^2/abc(1)^2+...
        (pcloud(:,2)-rc(2)*ovec).^2/abc(2)^2+...
        (pcloud(:,3)-rc(3)*ovec).^2/abc(3)^2 )<=1;

% take only the points that comprise an entire element
nppe=size(bndele,2); % number of points per element 
pts_in=find(indx_in);
B=zeros(size(bndele));
for ii=1:size(B,2)
    B(:,ii)=ismember(bndele(:,ii),pts_in);
end
ckele=sum(B,2)==nppe;

elein=find(ckele);
