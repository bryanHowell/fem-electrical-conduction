
%%

tic;
filename='headmesh_O2.mphtxt';
meshdata=read_mphtxt(filename);
tloadmesh=toc;

% COMSOL order is not the same as what is desired
meshdata.volele(:,5:10)=meshdata.volele(:,[5,6,8,7,10,9]);

%% assess node order in tetrahedral element

% a=meshdata.volele(1,:);
% p=meshdata.nodes(a,:);
% i1=[1,2,3,1];
% i2=[1,3,4,1];
% i3=[1,2,4,1];
% i4=[2,3,4,2];
% t=5:7;
% tt=8:10;
% 
% hold on;
% plot3(p(i1,1),p(i1,2),p(i1,3),'k.--','markersize',15);
% plot3(p(i2,1),p(i2,2),p(i2,3),'b--');
% plot3(p(i3,1),p(i3,2),p(i3,3),'r--');
% plot3(p(i4,1),p(i4,2),p(i4,3),'g--');
% plot3(p(5,1),p(5,2),p(5,3),'ks');
% plot3(p(6,1),p(6,2),p(6,3),'bs');
% plot3(p(7,1),p(7,2),p(7,3),'gs');
% plot3(p(8,1),p(8,2),p(8,3),'ko');
% plot3(p(9,1),p(9,2),p(9,3),'bo');
% plot3(p(10,1),p(10,2),p(10,3),'go');
% hold off;
% view([1,-1,1]);
% axis square;

%% find electrode contacts

p1=[79.3501;124.109;108.927]; % centroid of contact 1
bndele1=find_bndele(meshdata.nodes,meshdata.bndele,p1,[1;1;1]);
p2=[86.4755;131.8416;0]; % bottom of neck center
bndele2=find_bndele(meshdata.nodes,meshdata.bndele,p2,[30;30;5]);

atmp=meshdata.volele(:,1:3);
aa=unique(atmp(:));
btmp=meshdata.bndele(bndele1,:);
bb=unique(btmp);
ctmp=meshdata.bndele(bndele2,:);
cc=unique(ctmp);

hold on;
plot3(meshdata.nodes(aa,1),meshdata.nodes(aa,2),...
    meshdata.nodes(aa,3),'k.','markersize',2);
plot3(meshdata.nodes(bb,1),meshdata.nodes(bb,2),...
    meshdata.nodes(bb,3),'r.','markersize',15);
plot3(meshdata.nodes(cc,1),meshdata.nodes(cc,2),...
    meshdata.nodes(cc,3),'b.','markersize',15);
hold off;


