clc;

%% read data

xyzdata=load('dcfpath_p1_300axD3um.txt');
p1=[79.3501;124.109;108.927]; % centroid of contact 1

%% extract coordinates

numax=max(xyzdata(:,1));
xyz_ax=cell(numax,1);
xyz_node=cell(numax,1);

for k=1:numax
    ck=xyzdata(:,1)==k;
    xyz_ax{k}=xyzdata(ck,2:4)';
    npts=size(xyz_ax{k},2);
    xyz_node{k}=xyz_ax{k}(:,1:8:npts);
end

%% distance distributions

numnode=size(xyz_node{1},2);

r_ef=zeros(numax,1);
for k=1:numax
    rall=sqrt(sum((xyz_node{k}-p1*ones(1,numnode)).^2,1));
    r_ef(k)=min(rall);
end

[~,isort1]=sort(-r_ef(1:100));
[~,isort2]=sort(-r_ef(101:200));

%% split axons into two populations
% - 100 axons for CST
% - 100 axons for Hyperdirect pathway

xyz_cst=[];
xyz_hd=[];

for k=1:80
    
    kk=100+k;
    
    tmp1=[k*ones(numnode,1),xyz_node{isort1(k)}'];
    tmp2=[k*ones(numnode,1),xyz_node{100+isort2(k)}'];
    
    xyz_cst=[xyz_cst;tmp1];
    xyz_hd=[xyz_hd;tmp2];
    
end

hold on;
plot3(xyz_cst(:,2),xyz_cst(:,3),xyz_cst(:,4),'k.');
plot3(xyz_hd(:,2),xyz_hd(:,3),xyz_hd(:,4),'r.');
plot3(p1(1),p1(2),p1(3),'g.','markersize',20);
hold off;

dlmwrite('cst_p1_80axD3um.txt',xyz_cst);
dlmwrite('hd_p1_80axD3um.txt',xyz_hd);

%% visualize

figure;
hold on;
for k=1:numax
    if(k<=100)
        plot3(xyz_node{k}(1,:),xyz_node{k}(2,:),xyz_node{k}(3,:),'k.');
    elseif(k>100&&k<=200)
        plot3(xyz_node{k}(1,:),xyz_node{k}(2,:),xyz_node{k}(3,:),'r.');        
    else
%         plot3(xyz_node{k}(1,:),xyz_node{k}(2,:),xyz_node{k}(3,:),'r.');   
    end
end
hold off; 
