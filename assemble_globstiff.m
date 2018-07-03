function [A] = assemble_globstiff()

global meshdata; % mesh data
global K; % tensor field of material properties

%% assemble stiffness matrix

Ae=zeros(meshdata.ns,meshdata.ns,meshdata.nele); % preallocate

% coordinates (1st) for each shape function (2nd) for all elements (3rd)
coord=zeros(meshdata.nd,meshdata.ns,meshdata.nele);
for ii=1:meshdata.nd
    for jj=1:meshdata.ns
        coord(ii,jj,:)=meshdata.nodes(meshdata.volele(:,jj),ii);        
    end
end

for k=1:meshdata.nqp
    
    % shape function derivatives (1st) for each shape function (2nd)
    dshape=dirshape_ele(meshdata.gquad(:,k));

    % calculate Jacobian (1st and 2nd) for all elements (3rd)
    J=calc_jacob(dshape,coord);

    % calcuclate inverse and determinate of Jacobian
    [invJ,detJ]=calc_inv(J);

    % calculate shape function derivatives in physical domain
    gradN=calc_dirshape(dshape,invJ);

    % apply material properties to calculate flux
    % - map gradN to flux space
    KgradN=calc_flux(K,gradN);

    % assemble element stiffness matrix
    Ae=Ae+meshdata.gquad_weights(k)*(-1)*calc_elestiff(gradN,KgradN,detJ);

end

% assemble global stiffness matrix
Y=reshape(repmat(meshdata.volele,1,meshdata.ns)',...
    meshdata.ns,meshdata.ns,meshdata.nele);
X=permute(Y,[2,1,3]);
A=sparse(X(:),Y(:),Ae(:));

