clc;

%% load data

data1=load('phi_homiso.txt');
data2=load('phi_hetiso.txt');
data3=load('phi_hetani.txt');

%% visualize data

figure;
hold on;
plot(data1(:,3),data1(:,4),'k','linewidth',3);
plot(data2(:,3),data2(:,4),'b','linewidth',3);
plot(data3(:,3),data3(:,4),'r','linewidth',3);
hold off;
xlabel('z displacement (mm)','fontsize',30);
ylabel('potential (V)','fontsize',30);
set(gca,'fontsize',24);
axis tight;
legend('homogeneous, isotropic','heterogeneous, isotropic',...
    'heterogeneous, anisotropic','location','nw');
