function [A] = assemble_globstiff()

global meshdata; % mesh data
global K; % tensor field of material properties

%%

A=sparse(meshdata.ndof,meshdata.ndof);
Ae=zeros(meshdata.ns,meshdata.ns); % preallocate
Ke=zeros(meshdata.nd,meshdata.nd);

for e=1:meshdata.nele
    
    Ae=0*Ae;
    Ke=0*Ke;
    
    for ii=1:meshdata.nqp
        
        % current quadrature point
        qp_ii=meshdata.gquad(:,ii);
        
        % points making up the current element
        xyz=meshdata.nodes(meshdata.volele(e,:),:);
        
        % sfxn_data, shape function data
        % - numnodes by 4
        % - three columns: 
        % * c1 = shape function at ele. coordinates
        % * c2 = sfxn derivative w/r to psi ""
        % * c3 = sfxn derivative w/r to neu ""
        % * c4 = sfxn derivative w/r to zeta ""
        sfxn_data=shapefxn_3DO2(qp_ii(1),qp_ii(2),qp_ii(3));

        % derivatives of xy coord. w/r to ele. coord.
        dxde=sfxn_data(:,2)'*xyz(:,1); % scalar
        dxdn=sfxn_data(:,3)'*xyz(:,1); % ""
        dxdc=sfxn_data(:,4)'*xyz(:,1); % ""
        dyde=sfxn_data(:,2)'*xyz(:,2); % ""
        dydn=sfxn_data(:,3)'*xyz(:,2); % ""
        dydc=sfxn_data(:,4)'*xyz(:,2); % ""
        dzde=sfxn_data(:,2)'*xyz(:,3); % ""
        dzdn=sfxn_data(:,3)'*xyz(:,3); % ""
        dzdc=sfxn_data(:,4)'*xyz(:,3); % ""
        
        % Jacobian and its determinate and inverse
        J=[dxde,dyde,dzde;dxdn,dydn,dzdn;dxdc,dydc,dzdc]; % 3x3
        detJ=det(J); % scalar
        Jinv=inv_3by3(J); % 3x3
        
        gradN=Jinv*sfxn_data(:,2:4)'; % 3x6
        
        Ke=K(:,:,e);
%         Ke([1,4,5,7,8,9])=K(e,:); % upper triangle + diag.
%         Ke([2,3,6])=K(e,[2,4,5]); % lower triangle (symmetry)
                
        % element stiffness is 6x6
        Ae=Ae+meshdata.gquad_weights(ii)*(-1)*gradN'*(Ke*gradN)*detJ;
        
    end
    
%     li_A=nbyn_comb(meshdata.volele(e,:),numdof);
%     A(li_A)=A(li_A)+Ae(:);
    
    for a=1:meshdata.ns
        for b=1:meshdata.ns
            ag=meshdata.volele(e,a);
            bg=meshdata.volele(e,b);
            A(ag,bg)=A(ag,bg)+Ae(a,b);
        end
    end

end

