clear all
clc


load('data\receptor_data_glasser.csv')
load('data\SDI_conte69')
%receptor_names = readNPY('H:\科研\data_修改\neuromap_test\hansen_receptors-main\data\receptor_names_pet.npy');

data_g=ft_read_cifti('data\Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
% g1=data_g.gradient_1(1:64984);

posL=find(data_g.brainstructure==1);
posR=find(data_g.brainstructure==2);

SDI_conte69_L=SDI_conte69(posL);
SDI_conte69_R=SDI_conte69(posR);
[surf_lh, surf_rh] = load_conte69;
[sphere_lh, sphere_rh] = load_conte69('spheres');


for i=1:360
    pos=find(data_g.indexmax==i);
    for j = 1: size(receptor_data_glasser,2)
        receptor_data = receptor_data_glasser(:,j);
        receptor_data_32k(pos,j)=receptor_data(i);
    end
end




for i = 1 :size(receptor_data_32k,2)
    
    
    % load the conte69 hemisphere surfaces and spheres
    
    % Let's create some rotations
    n_permutations = 1000;
    y_rand = spin_permutations({SDI_conte69_L,SDI_conte69_R}, ...
        {sphere_lh,sphere_rh}, ...
        n_permutations,'random_state',0);
    
    SDI_rotated = squeeze([y_rand{1}(:,1,:); y_rand{2}(:,1,:)]);
    % for j = 1:n_permutations
    %     for i=1:360
    %     pos=find(gla360conte69==i);
    %     SDI_rotated_360(i,j)=mean(SDI_rotated(pos,j));
    %     end
    % end
    
    [r_original(i), pval(i)] = corr(receptor_data_32k(:,i),SDI_conte69,'rows','pairwise','type','spearman');
    r_rand(:,i) = corr(receptor_data_32k(:,i),SDI_rotated,'rows','pairwise','type','spearman');
    
    % [r_original_SDI, pval_SDI] = corr(g1_360,mean_SDI,'rows','pairwise','type','spearman');
    %  r_rand_SDI = corr(g1_360,SDI_rotated_360,'rows','pairwise','type','spearman');
    
    % Compute percentile rank.
    
    prctile_rank(i) = mean(r_original(i) > r_rand(:,i));
    significant(i) = prctile_rank(i) < 0.025 || prctile_rank(i) >= 0.975;
    
end


receptor_name = table2cell(readtable('data\receptor_name.xlsx'));
receptor_sig_name = receptor_name(find(significant == 1));
prctile_rank_sig = prctile_rank(find(significant == 1));
% gla360conte69 = data_g.indexmax;

r_sig = r_original(find(significant == 1));
r_p = [r_sig;prctile_rank_sig];

receptor_sig = receptor_data_glasser(:,(find(significant == 1)));

% for j=1:360
%     pos=find(gla360conte69==j);
%     receptor__sig_glasser(j,:) = mean(neuromap_sig(pos,:));
% end


%% Figure

% X1=sort(mean_SDI);
% Y1=g1_360;
load('data/mean_SDI.mat');

mean_SDI=log(mean_SDI);

X = mean_SDI;

% X = SDI_conte69;
% Y_receptor_sig = receptor_data_32k(:,[10,15,16,17]);
Y_receptor_sig = receptor_data_glasser(:,[10,15,16,17]);
nRe = size(Y_receptor_sig,2);

N = length(mean_SDI);
% N = length(SDI_conte69);%得到序列的长度
Xrank = zeros(N , 1); %存储X中各元素的排行
Yrank = zeros(N, 1); %存储Y中各元素的排行

%计算Xrank中的各个值
for i = 1 : N
    cont1 = 1; %记录大于特定元素的元素个数
    cont2 = -1; %记录与特定元素相同的元素个数
    for j = 1 : N
        if X(i) < X(j)
            cont1 = cont1 + 1;
        elseif X(i) == X(j)
            cont2 = cont2 + 1;
        end
    end
    Xrank(i) = cont1 + mean([0 : cont2]);
end

%计算Yrank中的各个值
for recp = 1:nRe
    Y = Y_receptor_sig(:,recp);
    for i = 1 : N
        cont1 = 1; %记录大于特定元素的元素个数
        cont2 = -1; %记录与特定元素相同的元素个数
        for j = 1 : N
            if Y(i) < Y(j)
                cont1 = cont1 + 1;
            elseif Y(i) == Y(j)
                cont2 = cont2 + 1;
            end
        end
        Yrank(i) = cont1 + mean([0 : cont2]);
    end
    
    
    X1 = Xrank;
    Y1 = Yrank;
    
    
    alpha = 0.01;
    [p1,s1]=polyfit(X1,Y1,1);
    UU1=linspace((min(X1)),(max(X1)));
    VV1=polyval(p1,UU1);
    [Y11,DELTA1] = polyconf(p1,UU1,s1,'alpha',alpha,'predopt','curve');
    
    r1=corrcoef(X1,Y1);
    r(recp)=r1(2,1);
    
    
    figure(recp);
    xlim([min(X1),max(X1)]);
    color = [0,0.5,0.5];
    h1=plot(X1,Y1,'ko');
    xlim([min(X1)-10,max(X1)]+10);ylim([min(X1)-10,max(X1)]+10);
    hold on
    h11=plot(UU1,Y11,UU1,Y11+DELTA1,'--',UU1,Y11-DELTA1,'--');
    hold on; set(h1(1),'MarkerFaceColor',[0,0,0],'MarkerSize',4,'Color',[1,1,1]);
    set(h11(1),'LineWidth',2,'Color',color);
    hold on;
    set(h11(2),'lineWidth',1,'Color',color);
    set(h11(3),'lineWidth',1,'Color',color);
    
    hold on;
    XX1=[UU1';flipud(UU1')];
    YY1=[(Y11+DELTA1)';flipud((Y11-DELTA1)')];
    fill(XX1,YY1,color,'facealpha',0.03,'edgealpha',0);
end

for recp = 1:nRe
    figure(recp+4);
    % histogram(r_rand_SDI,100);
    % set(gca,'xtick',[-0.3:0.1:0.3]);
    % hold on;
    bar(r_sig(recp),30,0.1/16);
    hold on;
    histogram(r_rand(:,recp),100);
    set(gca,'xtick',[-0.4:0.1:0.6]);
end

