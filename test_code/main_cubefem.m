clear all;
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

filename='cube_3dmesh_o2.vtk'; % O2 = order 2, quadratic
meshdata=readvtk_gmsh(filename);

% re-order points 5-10
% order of points from gmsh does not match the standard numberings as
% outlined in The Finite Element Method in Electromagnetics (by Jianming
% Jin)
% gmsh 5-10 -> 5,7,8,6,10,9
% COMSOL 5-10 -> 5,6,8,7,10,9
meshdata.volele(:,5:10)=meshdata.volele(:,[5,7,8,6,10,9]);
% meshdata.volele(:,5:10)=meshdata.volele(:,[5,6,8,7,10,9]);

%% define boundary conditions

surfnodes=unique(meshdata.bndele(:));
p_surf=meshdata.nodes(surfnodes,:);

% assign Dirichlet and Neumann BCs
ck_dir1=(p_surf(:,1)>=0.25&p_surf(:,1)<=0.75)...
        &(p_surf(:,2)>=0.25&p_surf(:,2)<=0.75)...
        &(p_surf(:,3)<0.5);
ck_dir2=(p_surf(:,1)>=0.25&p_surf(:,1)<=0.75)...
        &(p_surf(:,2)>=0.25&p_surf(:,2)<=0.75)...
        &(p_surf(:,3)>=0.5);
ck_neu=(~ck_dir1)&(~ck_dir2);

% % visualize BCs
% hold on;
% plot3(meshdata.nodes(surfnodes,1),meshdata.nodes(surfnodes,2),...
%       meshdata.nodes(surfnodes,3),'k.');
% plot3(meshdata.nodes(ck_dir1,1),meshdata.nodes(ck_dir1,2),...
%       meshdata.nodes(ck_dir1,3),'ro');
% plot3(meshdata.nodes(ck_dir2,1),meshdata.nodes(ck_dir2,2),...
%       meshdata.nodes(ck_dir2,3),'bs');
% plot3(meshdata.nodes(ck_neu,1),meshdata.nodes(ck_neu,2),...
%       meshdata.nodes(ck_neu,3),'go');  
% hold off;

meshdata.bndnodes=surfnodes;

% types of boundaries:
% Dirichlet: u=X, denoted by 1
% Neumann: -S*grad(u)=X, denoted by 2
meshdata.bndtype=zeros(length(surfnodes),1);
meshdata.bndtype(ck_dir1|ck_dir2)=1;
meshdata.bndtype(ck_neu)=2;

% ground is X=0 (dir1)
% insulation is -S*grad(u)=0 (neu)
% - 0 is default => no need to assign values for the above conditions
meshdata.bndval=zeros(length(surfnodes),1);
meshdata.bndval(ck_dir1)=1;

% clear unused variables
clear surfnodes;
clear p_surf;
clear ck_dir1;
clear ck_dir1;
clear ck_neu;

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

v1=([1,1,1]/norm([1,1,1]))';
v2=([-1,-1,2]/norm([-1,-1,2]))';
v3=cross(v1,v2);
V=[v1,v2,v3];
U=diag([0.9,0.1,0.1])*1e-3;
S=V*U*V';

K=S(:,:,ones(meshdata.nele,1));

% for ii=1:meshdata.ndof
%     K(ii,:)=S([1,4,5,7,8,9]);
% end

%% assemble system of equations

% Au=f
tic;
A=vecassembly_globstiff();
% force vector
% - no sources
f=zeros(size(meshdata.nodes,1),1);

tassemb=toc;

%% enforce Dirichlet boundary conditions

[A,f]=enforce_dirichletbc(A,f);

%% solve system

u=A\f;

scatter3(meshdata.nodes(:,1),meshdata.nodes(:,2),...
    meshdata.nodes(:,3),12,u,'filled');
xlabel('x (mm)','fontsize',30);
ylabel('y (mm)','fontsize',30);
zlabel('z (mm)','fontsize',30);
set(gca,'fontsize',24);
colorbar;
axis square;
colormap('jet');

return

%%

% aa=meshdata.nodes;
% isurf=unique(meshdata.bndele(:));
% 
% ea=[1,1,1,2,3,4];
% eb=[2,3,4,3,4,2];
% 
% % plot3(aa(:,1),aa(:,2),aa(:,3),'k.');
% % hold on;
% % plot3(aa(isurf,1),aa(isurf,2),aa(isurf,3),'g.');
% for ii=1:size(meshdata.intele,1)
%     
%     a=meshdata.nodes(meshdata.intele(ii,:),:); 
%     
%     for jj=5:10
%         
%         kk=jj-4;
%         
%         clf;
%         hold on;
%         plot3(a([1,2,3,1],1),a([1,2,3,1],2),a([1,2,3,1],3),'r-');
%         plot3(a([1,3,4,1],1),a([1,3,4,1],2),a([1,3,4,1],3),'r-');
%         plot3(a([2,3,4,2],1),a([2,3,4,2],2),a([2,3,4,2],3),'r-');
%         
%         plot3(a(1,1),a(1,2),a(1,3),'ko');
%         plot3(a(2,1),a(2,2),a(2,3),'bo');
%         plot3(a(3,1),a(3,2),a(3,3),'go');
%         plot3(a(4,1),a(4,2),a(4,3),'co');
%         
% %         plot3(a([ea(kk),eb(kk)],1),a([ea(kk),eb(kk)],2),...
% %             a([ea(kk),eb(kk)],3),'b-');
%         plot3(a(jj,1),a(jj,2),a(jj,3),'k.');
%         gtmp=input('');
%         hold off;
%         
%     end
%     
% end
% % hold off;

