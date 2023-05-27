clear all;
clc;

load('data/mean_SDI.mat');
mean_SDI = log2(mean_SDI);
% load('mean_SDI_thr.mat');
label = xlsread('data/glasser360_7networks.xlsx');
label = label(:,3);
mean_SDI_label=[label,mean_SDI];
% mean_SDI_thr(mean_SDI_thr==1)=0;
% mean_SDI_thr_label=[label,mean_SDI_thr];
mean_SDI_label_sort=sortrows(mean_SDI_label,1);

%%
y=xlsread('data/mean_SDI_7labels.xlsx');
figure;
boxplot(y,'colors','k','orientation', 'horizontal');
set(gca,'yticklabel',{'VIS','SM','DA','VA','LIM','FP','DM'});

for i = 1:7
    p=find(label==i);
    p_SDI=mean_SDI(p,:);
    for j = 1:7
        if i ~= j
            q=find(label==j);
            q_SDI=mean_SDI(q,:);
            [h(i,j),p_h(i,j)]=ttest2(p_SDI,q_SDI);
        end   
    end
    FDR(i,:) = mafdr(p_h(i,:),'BHFDR', true);
end

U=triu (FDR);
L = tril(FDR);
L(find(L==0))=1;

figure;
imagesc(abs(log2(L)));
colormap(othercolor('BuPu8'));
grid on;
% grid minor;
% hold on;
set(gca,'xticklabel',[],'yticklabel',[],'xtick',(0.5:1:7.5),'ytick',(0.5:1:7.5))

xtick = (0.5:1:7.5);
ytick = (0.5:1:7.5);




