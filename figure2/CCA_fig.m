clear all;
clc;



beha_coeff = xlsread('data\2beha_without_adjust.xlsx');
SDI_coeff = xlsread('data\2SDI_age_sex_FD.xlsx');

beha_coeff1 = beha_coeff(:,1);
SDI_coeff1 = SDI_coeff(:,1);


label = xlsread('data\glasser360_7networks.xlsx');
SDI_label = label(:,3);
beha_label=xlsread('data\beha_label.xlsx');
beha_label = beha_label(:,1);
beha_label = beha_label((1:108),:);

beha_cl_1 = [beha_coeff1,beha_label];
SDI_cl_1 = [SDI_coeff1,SDI_label];
SDI_cl_1_r = sortrows(SDI_cl_1,2);
%% mode 1
b1= SDI_cl_1_r(:,1);
a1 = beha_cl_1(:,1);
x1=1:1:size(b1);
x2=1:1:size(a1);
figure;


scatter(x1(1:59),b1(1:59),'r','filled','MarkerEdgeColor',[0.5 0 0],'LineWidth',1.5);
hold on;
scatter(x1(60:112),b1(60:112),'b','filled','MarkerEdgeColor',[0 0 .5],'LineWidth',1.5);
hold on;
scatter(x1(113:156),b1(113:156),'y','filled','MarkerEdgeColor',[0.5 0.5 0],'LineWidth',1.5);
hold on;
scatter(x1(157:204),b1(157:204),'g','filled','MarkerEdgeColor',[0 .5 0],'LineWidth',1.5);
hold on;
scatter(x1(205:233),b1(205:233),'c','filled','MarkerEdgeColor',[0 .5 .5],'LineWidth',1.5);
hold on;
scatter(x1(234:278),b1(234:278),'MarkerFaceColor',[0.4 0.3 0.8],'MarkerEdgeColor',[0.2 0.2 0.5],'LineWidth',1.5);
hold on;
scatter(x1(279:360),b1(279:360),'m','filled','MarkerEdgeColor',[0.5 0 .5],'LineWidth',1.5);

figure;
scatter(x2(1:6),a1(1:6),'r','filled','MarkerEdgeColor',[0.5 0 0],'LineWidth',1.5);
hold on;
scatter(x2(7:30),a1(7:30),'b','filled','MarkerEdgeColor',[0 0 .5],'LineWidth',1.5);
hold on;
scatter(x2(31:54),a1(31:54),'y','filled','MarkerEdgeColor',[0.5 0.5 0],'LineWidth',1.5);
hold on;
scatter(x2(55:58),a1(55:58),'g','filled','MarkerEdgeColor',[0 .5 0],'LineWidth',1.5);
hold on;
scatter(x2(59:63),a1(59:63),'c','filled','MarkerEdgeColor',[0 .5 .5],'LineWidth',1.5);
hold on;
scatter(x2(64:67),a1(64:67),'MarkerFaceColor',[0.4 0.3 0.8],'MarkerEdgeColor',[0.2 0.2 0.5],'LineWidth',1.5);
hold on;
scatter(x2(68:85),a1(68:85),'m','filled','MarkerEdgeColor',[0.5 0 .5],'LineWidth',1.5);
hold on;
scatter(x2(86:102),a1(86:102),'MarkerFaceColor',[0.7 0.7 0],'MarkerEdgeColor',[0.5 0.5 0],'LineWidth',1.5);
hold on;
scatter(x2(103:108),a1(103:108),'MarkerFaceColor',[0.9 0.6 0.3],'MarkerEdgeColor',[0.8 0.4 0.1],'LineWidth',1.5);






