function [gradN_phy] = calc_dirshape(gradN_ele,invJ)
% gradN_ele, derivatives of shape functions in parent domain (nd x ns)
% invJ, inverse of Jacobian for all elements (nd x nd x ne)
% gradN_phy,derivatives of shape functions in physical domain (nd x ns x ne) 
%
% nd = # of dimensions
% ns = # of shape functions
% ne = # of elements

[nd,~,ne]=size(invJ);
ns=size(gradN_ele,2);

gradN_phy=zeros(nd,ns,ne); % pre-allocate

for jj=1:ns
    
    gNejj=gradN_ele(:,jj)';
    N=gNejj(ones(nd,1),:,ones(ne,1));
    gradN_phy(:,jj,:)=sum(invJ.*N,2);    
    
end