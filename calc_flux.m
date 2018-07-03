function [flux_u] = calc_flux(K,grad_u)
% K, tensor field of material properties (nd x nd x ne)
% grad_u, gradient of quantity (nd x ns x ne)
% flux_u, flux of quantity (nd x ns x ne)
%
% nd = # of dimensions
% ns = # of shape functions
% ne = # of elements

[nd,ns,ne]=size(grad_u);
flux_u=zeros(nd,ns,ne); % pre-allocate

for jj=1:ns
    
    gu_jj=reshape(grad_u(:,jj,:),1,nd,ne);
    U=gu_jj(ones(nd,1),:,:);
    flux_u(:,jj,:)=sum(K.*U,2);
        
end