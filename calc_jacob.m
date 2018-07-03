function [J] = calc_jacob(gradN,P)
% gradN, derivative of shape functions (nd x ns)
% P, collection of points per element (nd x ns x ne)
% J, Jacobian of element in parent (i.e., element) domain (nd x nd x ne)
%
% nd = # of dimensions
% ns = # of shape functions
% ne = # of elements

[nd,~,ne]=size(P);

J=zeros(nd,nd,ne); % pre-allocate

for ii=1:nd
    
    dNdii=gradN(ii,:);
    N=dNdii;
    N=N(ones(nd,1),:,ones(ne,1));    
    
    dxyz_dii=sum(N.*P,2);
    J(ii,:,:)=reshape(dxyz_dii,1,nd,ne);
    
end
