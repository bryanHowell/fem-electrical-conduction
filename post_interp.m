function [u_int] = post_interp(pint,u_fem)
% u_fem, FEM solution at all degrees of freedom (i.e., nodes) 
% pint, set of points to interpolate
% u_int, solution at interpolated points

% short description:
% 1. find the parent element
% - the point to be interpolated will reside within one element
% - locate which element contains the point
% 2. inverse mapping
% - we have physical coordinates but the FEM uses parent coordinates
% - physical coordinates = xyz and parent coordinates = enc
% - we need to find which enc coordinate maps to xyz (an inverse problem)
% 3. interpolate
% - use shape functions to interpolate the solution in parent space
%
% nd = # of dimensions
% ns = # of shape functions
% ne = # of elements

global meshdata;

%% pre-processing
% there may be 100s of thousands of elements or more
% (research-grade problems typically require millions of elements)
% => need to reduce the number of elements being searched

% first, bookkeeping...
% elecoord, holds all points for all shape fxns for all elements
elecoord=zeros(meshdata.nd,meshdata.ns,meshdata.nele);
for ii=1:meshdata.nd
    for jj=1:meshdata.ns
        elecoord(ii,jj,:)=meshdata.nodes(meshdata.volele(:,jj),ii);        
    end
end

% rectangular prisms that circumscribe all elements
% nd x 1 x ne -> squeeze -> nd x ne
% (squeeze removes singular dimensions)
cbox=[min(elecoord,[],2),max(elecoord,[],2)];

% define rectangular prism that bounds interpolated points
xyz_min=min(pint,[],2);
xyz_max=max(pint,[],2);
% the user may search for a single point, points along a line, or points in
% in a place. => some dimensions may be singular
cksing=(xyz_max-xyz_min)==0;
if(any(cksing))
    maxdel=max(cbox(:,2,:)-cbox(:,1,:),[],3);
    % 11/20*a - -11/20*a = 11/10*a => safety factor of 1.1
    xyz_min(cksing)=xyz_min(cksing)-(11/20)*maxdel(cksing);
    xyz_max(cksing)=xyz_max(cksing)+(11/20)*maxdel(cksing);
end

% check which elements are within the rect. prism that bounds the interp.
% points
ckele=zeros(meshdata.nd,4,meshdata.nele);
for ii=1:meshdata.nd
    % only the first four columns need to be used as these vertices define
    % the tetrahedron. The rest of the points are on the edges of the tet.
    ckele(ii,:,:)=elecoord(ii,1:4,:)>=xyz_min(ii)...
        &elecoord(ii,1:4,:)<=xyz_max(ii);
end
ck_xyzin=sum(ckele,1)==3; % if xyz in prism, sum is 3 across rows
ck_elein=sum(ck_xyzin,2)>0; % if any points in ele. are w/in prism, sum>0
ck_elein=reshape(ck_elein,meshdata.nele,1,1);

subele=find(ck_elein);
num_subele=length(subele);

%% stage 1: find parent element

numpts=size(pint,2);
pntele=zeros(numpts,1);

cboxsub=cbox(:,:,subele);

for k=1:numpts
    
    p=pint(:,k);
    pp=p(:,:,ones(num_subele,1)); % repeat points in 3rd dimension
    
    cktmp=pp>=cboxsub(:,1,:)&pp<=cboxsub(:,2,:); % points in boxes?
    ckpts=reshape(sum(cktmp,1)==3,num_subele,1,1);

    candele=subele(ckpts); % elements that are candidates

    for ii=1:length(candele)
        
        vert=elecoord(:,1:4,candele(ii));
        ck=insidetet(vert(:,1),vert(:,2),vert(:,3),vert(:,4),p);
        if(ck==1)
            pntele(k)=candele(ii);
            break;            
        end    
        
    end
    
    % if no candidate successful, look for nearest neighbour amongst
    % candidates
    if(pntele(k)==0)
        
        tic;
        numcand=length(candele);
        rcand=squeeze( mean(elecoord(:,1:4,candele),2) );
        rall=sqrt(sum((rcand-p*ones(1,numcand)).^2,1));
        [~,inn]=min(rall);
        pntele(k)=candele(inn);
        
    end

end

%% stage 2: inverse mapping

pntcoord=zeros(meshdata.nd,numpts);
epstol=1e-8;
maxcnt=10;
count=0;
for k=1:numpts
    
    p=pint(:,k); % current point (nd x 1)

    elexyz=elecoord(:,:,pntele(k)); % points comprising ele. (nd x ns)
    
    % initial guess
    if(k==1)
        % neutral point (0<=e,n,c<=1) =>...
        enc=(1/2)*ones(meshdata.nd,1);
    else
        % next point is likely to be by previous point
        enc=pntcoord(:,k-1);
    end
    
    R=1e3;
    count=0;
    while(norm(R)>epstol&&count<maxcnt)
        
        N=shapefxn_tet(enc); % shape fxn @ (e,n,c) (1 x ns)
        gradN=dirshape_ele(enc); % grad of shape fxn @ (e,n,c) (nd x ns)

        R=p-(N*elexyz')'; % residual vector (nd x 1)
        J=-(gradN*elexyz')'; % Jacobian (nd x nd)
        Jinv=calc_inv(J); % inverse of Jacobian (nd x nd)

        enc=enc-Jinv*R; % update guess     
        count=count+1;
        if(count>=maxcnt)
            display('welp!');
        end
    end

    pntcoord(:,k)=enc;
    
end

%% stage 3: interpolate

u_int=zeros(1,numpts);

for k=1:numpts
    
    N=shapefxn_tet(pntcoord(:,k)); % 1 x ns
    uele=u_fem(meshdata.volele(pntele(k),:)); % ns x 1
    u_int(k)=N*uele; % interpolated point in parent element
    
end

