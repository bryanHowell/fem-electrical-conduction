function [Aele] = calc_elestiff(gradN,fluxN,detJ)
% gradN, gradient of shape functions (nd x ns x ne)
% fluxN, material tensor * gradN (nd x ns x ne)
% detJ, determinate of jacobian (1 x ne)
%
% nd = # of dimensions
% ns = # of shape functions
% ne = # of elements

[~,ns,ne]=size(gradN);
Aele=zeros(ns,ns,ne);% pre-allocate

for ii=1:ns
    
    fNjj=fluxN(:,ii,:);
    F=fNjj(:,ones(ns,1),:);
    Aele(ii,:,:)=sum(gradN.*F,1);
        
end

DJ=reshape(detJ,1,1,ne);
DJ=DJ(ones(ns,1),ones(ns,1),:);
Aele=Aele.*DJ;

