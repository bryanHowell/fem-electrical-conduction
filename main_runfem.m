% clear all;
clc;

%% global variables
% - mesh data is very large
% - MATLAB cannot pass values by reference
% - passing values via a copy will be slow
% => although not good practice, in general, I will a few global
% declarations to circumvent the above

global meshdata;
global K;

%% load mesh
% - mesh was created using gmsh

filename='headmesh_O2.mphtxt'; % O2 = order 2, quadratic
tic;
meshdata=read_mphtxt(filename);
tmesh=toc;

% re-order points 5-10
% order of points from gmsh does not match the standard numberings as
% outlined in The Finite Element Method in Electromagnetics (by Jianming
% Jin)
% gmsh 5-10 -> 5,7,8,6,10,9
% COMSOL 5-10 -> 5,6,8,7,10,9
% meshdata.volele(:,5:10)=meshdata.volele(:,[5,7,8,6,10,9]);
meshdata.volele(:,5:10)=meshdata.volele(:,[5,6,8,7,10,9]);

display(['mesh time: ',num2str(tmesh),' s']);

%% define boundary conditions

p1=[79.3501;124.109;108.927]; % centroid of contact 1
bndele1=find_bndele(meshdata.nodes,meshdata.bndele,p1,[1;1;1]);
p2=[86.4755;131.8416;0]; % bottom of neck center
bndele2=find_bndele(meshdata.nodes,meshdata.bndele,p2,[30;30;5]);

meshdata.setbndele=[bndele1;bndele2];
meshdata.setbndele_type=ones(size([bndele1;bndele2])); % 1 = Dirichlet
meshdata.setbndele_val=[ones(size(bndele1));zeros(size(bndele2))];

%% Guassian quadrature 
% - parameters for numerical integration
a=(5+3*sqrt(5))/20;
b=(5-sqrt(5))/20;
meshdata.gquad=[a,b,b,b;b,a,b,b;b,b,a,b];
meshdata.gquad_weights=[1,1,1,1]/24;

%% define counts

meshdata.nd=3; % # of dimensions
meshdata.ndof=size(meshdata.nodes,1); % # of degrees of freedom
meshdata.nele=size(meshdata.volele,1); % # of domain elements
meshdata.nbdn=size(meshdata.bndele,1); % # of boundary elements 
meshdata.ns=size(meshdata.volele,2); % # of shape functions
meshdata.nqp=size(meshdata.gquad,2); % # of quadrature points

%% define tensor field
% in 3D, K is 3x3
% K=[kxx,kxy,kxz;kyx,kyy,kyz;kzx,kzy,kzz]
% - K is symmetric (kxy=kyx, kxz=kzx, kyz=kzy)
% => only six values are needed

tic;

% tensor field
K=zeros(3,3,meshdata.nele);

% types of volume elements (see meshdata.volele_type)
% regions: 1=soft tissue, 2=skull, 3-4=brain, 5-6=glial scar
sigval=[0.33,0.02,0.23,0.23,0.13,0.13]*1e-3; % S/m -> S/mm

% HETEROGENEOUS AND ANISOTROPIC
% - anisotropic brain
% - isotropic skull and soft tissues

% read in anisotropic tensor field in brain
% organization of tfield
% - 9 columns
% - 1st 3 columns: xyz
% - last 6 columns: a11, a12, a22, a13, a23, a33
% (the above correspond to linear indices 1,4,5,7,8,9 in a 3x3 in MATLAB)

tfield=load('tfield_lpres_p1.txt');
n_tf=size(tfield,1);

ckbrain=meshdata.volele_type==3|meshdata.volele_type==4;
ibrain=find(ckbrain);
clear ckbrain;

% determine bounds of tensor field
xyztf_min=min(tfield(:,1:3),[],1);
xyztf_max=max(tfield(:,1:3),[],1);
[~,idim]=max(xyztf_max-xyztf_min);
% breaking dimension that spans greatest distance into smaller quadrants
% - this is essential, otherwise mapping data from tensor field to mesh
% takes on the order of hours!!!
nquad=100;
do=xyztf_min(idim);
df=xyztf_max(idim);
deld_bnds=do:(df-do)/nquad:df;
quad_pts=cell(nquad,1);
quad_ele=cell(nquad,1);
for k=1:nquad
    cksub=tfield(:,idim)>=deld_bnds(k)&tfield(:,idim)<=deld_bnds(k+1);
    quad_pts{k}=tfield(cksub,1:3)';
    quad_ele{k}=find(cksub);
end

% nearest neighbour approach
for k=1:length(ibrain)
      
    % current element number    
    kk=ibrain(k);
    
    % centroid of tetrahedron
    rcent=mean(meshdata.nodes(meshdata.volele(kk,1:4),:),1)';
    
    % find nearest neighbour in a given quadrant
    iquad=sum(rcent(idim)>=deld_bnds);
    nele=size(quad_pts{iquad},2);
    all_r=sqrt( sum((quad_pts{iquad}-rcent*ones(1,nele)).^2,1) ); 
    [~,imin]=min(all_r);
    inn=quad_ele{iquad}(imin); % nearest neighbour in tensor field
    
    S=zeros(3,3);
    S([1,4,5,7,8,9])=tfield(inn,4:9); % upper triangular + diagonal
    S([2,3,6])=S([4,7,8]); % aij = aji (symmetry)
    K(:,:,kk)=S*1e-3; %S/m -> S/mm
    
end

% define isotropic tensors for skull and brain
for ii=[1,2,5,6]
    
    ckreg=meshdata.volele_type==ii;
    volele_ii=find(ckreg);
    nii=length(volele_ii);
    
    S=diag(sigval(ii)*[1;1;1]);
    S=S(:,:,ones(nii,1));
    K(:,:,volele_ii)=S;
    
end

% % HETEROGENEOUS AND ISOTROPIC
% % - isotropic brain, skull, and soft tissues
% for ii=1:6
%     
%     ckreg=meshdata.volele_type==ii;
%     volele_ii=find(ckreg);
%     nii=length(volele_ii);
%     
%     S=diag(sigval(ii)*[1;1;1]);
%     S=S(:,:,ones(nii,1));
%     K(:,:,volele_ii)=S;
%     
% end

% % HOMOGENEOUS AND ISOTROPIC
% % - whole head treated like grey matter
%    
% S=diag(sigval(3)*[1;1;1]);
% S=S(:,:,ones(meshdata.nele,1));
% K=S;
%     
% ttensor=toc;
% display(['construct tensor field: ',num2str(ttensor),' s']);

%% assemble system of equations
% Au=f

tic;
A=assemble_globstiff();
% force vector
% - no sources
f=zeros(meshdata.ndof,1);
tassemb=toc;

display(['assembly time: ',num2str(tassemb),' s']);

%% enforce Dirichlet boundary conditions

tic;
[A,f]=enforce_dirichletbc(A,f);
tbc=toc;

display(['enforce boundary conditions: ',num2str(tbc),' s']);

%% solve system

tic;
u=A\f;
tsolve=toc;

display(['solve time: ',num2str(tsolve),' s']);

%% read axon coordinates

cd('model_analysis/');
cstdata=load('cst_p1_90axD3um.txt');
cd ..;

numax=max(cstdata(:,1));
xyz_cst=cell(numax,1);
for k=1:numax
    ck=cstdata(:,1)==k;
    xyz_cst{k}=cstdata(ck,2:4)';
end
numpts=size(xyz_cst{1},2);

xyzint=[];
for k=1:numax
    xyzint=[xyzint,xyz_cst{k}];
end

%% interpolation

% pc=mean(meshdata.nodes)';
% % xa=-50:0.01:50;
% % ya=1.5*ones(size(xa));
% % za=zeros(size(xa));
% za=-75:0.1:50;
% xa=3*ones(size(za));
% ya=3*ones(size(za));
% xyzint=[xa;ya;za]+p1*ones(1,length(xa));

tic;
u_int=post_interp(xyzint,u);
tinterp=toc;

U_int=zeros(numpts,numax);
U_int(:)=u_int;
U_int=U_int';

return

hold on;
plot3(meshdata.nodes(:,1),meshdata.nodes(:,2),meshdata.nodes(:,3),'k.',...
    'markersize',2);
plot3(xyzint(1,:),xyzint(2,:),xyzint(3,:),'b.','linewidth',3);
plot3(p1(1),p1(2),p1(3),'r.','markersize',15);
hold off;

%%

sub=u>0.01;

scatter3(meshdata.nodes(sub,1),meshdata.nodes(sub,2),...
    meshdata.nodes(sub,3),2,u(sub),'filled');
xlabel('x (mm)','fontsize',30);
ylabel('y (mm)','fontsize',30);
zlabel('z (mm)','fontsize',30);
set(gca,'fontsize',24);
colorbar;
axis square;
colormap('jet');
view([1,1,1]);
axis tight;

return

%%

aa=meshdata.points;
isurf=unique(meshdata.surfele(:));

ea=[1,1,1,2,3,4];
eb=[2,3,4,3,4,2];

% plot3(aa(:,1),aa(:,2),aa(:,3),'k.');
% hold on;
% plot3(aa(isurf,1),aa(isurf,2),aa(isurf,3),'g.');
for ii=1:size(meshdata.intele,1)
    
    a=meshdata.points(meshdata.intele(ii,:),:); 
    
    for jj=5:10
        
        kk=jj-4;
        
        clf;
        hold on;
        plot3(a([1,2,3,1],1),a([1,2,3,1],2),a([1,2,3,1],3),'r-');
        plot3(a([1,3,4,1],1),a([1,3,4,1],2),a([1,3,4,1],3),'r-');
        plot3(a([2,3,4,2],1),a([2,3,4,2],2),a([2,3,4,2],3),'r-');
        
        plot3(a(1,1),a(1,2),a(1,3),'ko');
        plot3(a(2,1),a(2,2),a(2,3),'bo');
        plot3(a(3,1),a(3,2),a(3,3),'go');
        plot3(a(4,1),a(4,2),a(4,3),'co');
        
%         plot3(a([ea(kk),eb(kk)],1),a([ea(kk),eb(kk)],2),...
%             a([ea(kk),eb(kk)],3),'b-');
        plot3(a(jj,1),a(jj,2),a(jj,3),'k.');
        gtmp=input('');
        hold off;
        
    end
    
end
% hold off;

