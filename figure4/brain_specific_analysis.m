
brain_specific_gene_id = readtable('data/ENSG2NCBIgeneID.txt');
brain_specific_gene_id = table2array(brain_specific_gene_id(:,2));

load('data/probeInformation.mat');
load('data/parcelExpression');
AHBA_EntrezID = probeInformation.EntrezID;
gene_name = probeInformation.GeneSymbol;

[c,ia,ib] = intersect(brain_specific_gene_id,AHBA_EntrezID);

gene = parcelExpression(:,(2:size(parcelExpression,2)));
gene_all = gene(all(~isnan(gene),2),:);
gene_specific = gene_all(:,ib);
gene_name = gene_name(ib,:);

% xlswrite('gene_name_brain_specific',gene_name,'sheet1');

load('data/SDI')
X = zscore(gene_specific);
Y = mean(SDI')';
Y = Y(181:360);
Y = zscore(Y(all(~isnan(gene),2)));

%% plsr

num = 15;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]= plsregress(X,Y,num);

% figure(2);
% plot(1:num,100*PCTVAR(2,:),'-bo','MarkerFaceColor','r');
% xlabel('Number of PLS components');
% ylabel('Expained variance');
% set(gcf,'unit','centimeters','position',[10 5 10 8]);

% 
X1=zscore(XS(:,1));
Y1=zscore(YS(:,1));
[R1,p1]=corr(X1,Y1);

%align PLS components with desired direction for interpretability
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end

PLS_Weights_W1 = stats.W(:,1);
[PLS1w,Indx1] = sort(PLS_Weights_W1,'descend');
XS1 = XS(:,1);
X_W = XL


%% Permutation testing
temp1=cumsum(100*PCTVAR(2,1:num));

ncomp=1;
Rsquared = temp1;
nIterNull=1000;
array_null_Rsq = zeros(nIterNull,1);
nnodes = size(X,1);
nterms = 1;
ngenes = size(X,2);


    for j=1:nIterNull
        Posp = randperm(nnodes)';
        Yp=Y(Posp);
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,num);
        temp=cumsum(100*PCTVAR(2,1:num));
        array_null_Rsq(j) = temp(ncomp);
        fprintf('第%d次置换检验\n',j);
    end


PLSR_pval = length(find(array_null_Rsq>=Rsquared))/nIterNull ;
PLSR_Rsq = Rsquared;

%% bootstrap

PLS1weights=[];
nBootstraps = 1000;
disp('... Bootstrapping ...')

for iB = 1 :nBootstraps
    
    fprintf('第%d次bootstrap\n',iB);
    myresample = randsample(size(X,1),size(X,1),1);
    % Resampling X,Y
    Xb = X(myresample,:);
    Yb = Y(myresample,:);
    
    [XLb,YLb,XSb,YSb,BETAb,PCTVARb,MSEb,statsb]= plsregress(Xb,Yb,num);
    tem=statsb.W(:,1);
    %     order the newly obtained weights the same way as initial PLS
    newWei=tem(Indx1);
    %     As the sign of PLS components is arbitrary, make sure this aligns between runs
    if corr(PLS1w,newWei)<0
        newWei=-1*newWei;
    end
    PLS1weights=[PLS1weights,newWei];
    
end

%standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');

%z-score weights
temp1=PLS1w./PLS1sw';
temp1_pos = find(abs(temp1)>2.58);


temp1_wei_boots = temp1(temp1_pos);
gene_name_boots = gene_name(temp1_pos);

[pos_boots_val,pos_boots] = sort(temp1_wei_boots(find(temp1_wei_boots>0)),'descend');
[neg_boots_val,neg_boots] = sort(temp1_wei_boots(find(temp1_wei_boots<0)),'descend');

[boots_val,boots] = sort(temp1_wei_boots,'descend');
[boots_val_abs,boots_abs] = sort(abs(temp1_wei_boots),'descend');

gene_name_sort_pos = gene_name_boots(pos_boots,:);
gene_name_sort_neg = gene_name_boots(neg_boots,:);
gene_name_sort = gene_name_boots(boots);
gene_name_abs_sort = gene_name_boots(boots_abs);

% xlswrite('gene_name_sort_pos',gene_name_sort_pos,'sheet1','A1:A554');
% xlswrite('gene_name_sort_pos',pos_boots_val,'sheet1','B1:B554');
% xlswrite('gene_name_sort_neg',gene_name_sort_neg,'sheet1','A1:A386');
% xlswrite('gene_name_sort_neg',neg_boots_val,'B1:B386');
% 
% xlswrite('gene_name_sort',gene_name_sort,'sheet1','A1:A940');
% xlswrite('gene_name_sort',boots_val,'sheet1','B1:B940');
% 
% xlswrite('gene_name_abs_sort',gene_name_abs_sort,'sheet1','A1:A940');
% xlswrite('gene_name_abs_sort',boots_val_abs,'sheet1','B1:B940');
% 
% xlswrite('brain_specific_gene_name',gene_name,'sheet1');

% alpha = 0.05;
% U11 = [ ones(length(Y),1), X1];
% [b,bint,r,rint,stats] = regress(Y, U11);
% p=stats(3);
% [p1,s1]=polyfit(X1,Y,1);
% UU1=linspace((min(X1)),(max(Y)));
% % xx1=linspace(-0.05,0.3);
% VV1=polyval(p1,UU1);
% [Y11,DELTA1] = polyconf(p1,UU1,s1,'alpha',alpha,'predopt','curve');
% 
% figure(3);
% color = [0,0.5,0.5];
% h1=plot(X1,Y,'ko');
% hold on
% h11=plot(UU1,Y11,UU1,Y11+DELTA1,'--',UU1,Y11-DELTA1,'--');
% hold on; set(h1(1),'MarkerFaceColor',[0,0,0],'MarkerSize',4,'Color',[1,1,1]);
% set(h11(1),'LineWidth',2,'Color',color);
% hold on;
% set(h11(2),'lineWidth',1,'Color',color);
% set(h11(3),'lineWidth',1,'Color',color);
% 
% hold on;
% XX1=[UU1';flipud(UU1')];
% YY1=[(Y11+DELTA1)';flipud((Y11-DELTA1)')];
% fill(XX1,YY1,color,'facealpha',0.03,'edgealpha',0);
% axis([(min(X1)-0.1) (max(X1)+0.1) -2.5 3]);
% set(gcf,'unit','centimeters','position',[10 5 10 8]);
% 
% for j=1:nIterNull
%     Posp = randperm(nnodes)';
%     Yp=Y(Posp);
%     r_p(j) = corr(X1,Yp);
%     fprintf('第%d次置换检验\n',j);
% end
% 
% pval = length(find(r_p>=R1))/nIterNull ;
% 
% figure;
% histogram(r_p,50)
% hold on;
% bar(R1,58,0.01);
% set(gca,'XLim',[-0.3,0.6]);


% permutation

%     cou(k) = sum(array_null_Rsq<=x_per(3)&array_null_Rsq>(3-1));
% FIGURE 4
% x_per = linspace(0,50);
% for k = 2:100
%     cou(k) = sum(array_null_Rsq<=x_per(k)&array_null_Rsq>x_per(k-1));
%     
% end


% figure(4);
% bar(cou);
% % set(gca,'XLim',[0,50]);
% 
% set(gca,'XTicklabel',(0:10:100))
% 
% hold on;
% bar(Rsquared*2,100);
% set(gcf,'unit','centimeters','position',[10 5 10 8]);







