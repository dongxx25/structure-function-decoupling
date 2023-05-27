%% split train set and test set

ID = readtable('data/family_ID.xlsx');

load('data/id_use.mat')
ID_all = table2array(ID(:,1));
[c, ia, ib] = intersect(ID_all, id_use);

Family_ID = table2cell(ID(:,2));
Family_ID = Family_ID(ia);


nSubj=size(Family_ID,1);
train_set=cell(10,10);
test_set=cell(10,10);
Famtrain_set=cell(10,10);
Famtest_set=cell(10,10);
cross=1;
uniqFam=unique(Family_ID);
nFam=size(uniqFam,1);
while cross<=10
    cross;
    kfoldout=10; % 10 CV for outer loop
    c_out = cvpartition(nFam,'Kfold',kfoldout); % for outer CV
    for fdout=1:10
        Famtrain_set{cross,fdout}=training(c_out,fdout);
        Famtest_set{cross,fdout}=test(c_out,fdout);
        temp=zeros(nSubj,1);
        temp(find(ismember(Family_ID,uniqFam(training(c_out,fdout)))))=1;
        train_set{cross,fdout}=logical(temp);
        temp=zeros(nSubj,1);
        temp(find(ismember(Family_ID,uniqFam(test(c_out,fdout)))))=1;
        test_set{cross,fdout}=logical(temp);
    end
    cross=cross+1;
end

%% predict behavioral scores
nSubj=877;nROI=360;

load('data/behavioral_data')

behav=behavioral_data(:,[13,14,27,28,30]);%cog

nItem=size(behav,2);
precorr=zeros(nItem,100);

label = xlsread('data/glasser360_7networks.xlsx');
label = label(:,3);
load('data/SDI');

vis = SDI(find(label==1),:);
sm = SDI(find(label==2),:);
da = SDI(label==3,:);
va = SDI(find(label==4),:);
lim = SDI(find(label==5),:);
fp = SDI(find(label==6),:);
dmn = SDI(find(label==7),:);


data=SDI;
% data=[da;va;lim;fp;dmn];
% data=[vis;sm];

%for
for out=1 :nItem
    %——————————————
    lab=behav(:,out)';
    
    acrr=zeros(10,10);
    % kfoldin=3;
    for cross=1%:10
        % kfoldout=10; % 10 CV for outer loop
        for fdout=1%:10
            fprintf('行为数据%d，第%d次10折交叉验证%d折\n',out, cross,fdout);
            Train_dat = data(:,train_set{cross,fdout});
            Train_lab = lab(train_set{cross,fdout});
            Test_dat = data(:,test_set{cross,fdout});
            Test_lab = lab(test_set{cross,fdout});
            
            
            % Test based on the optimal parameters obtained from the training data.
            
            %             indpout=1:size(Train_dat,1);
            %             [~,p_tmpout] = corr(Train_dat', Train_lab_new');
            %             indexout = find(p_tmpout<0.05);
            %             Train_dat = Train_dat(indexout,:);
            %             Test_dat = Test_dat(indexout,:);
            
            [B3{cross,fdout},FitInfo3{cross,fdout}]=lasso(Train_dat',Train_lab','Alpha',0.01,'CV',3);
            Train_lab_new = Train_dat'*B3{cross,fdout}(:,FitInfo3{cross,fdout}.IndexMinMSE) + FitInfo3{cross,fdout}.Intercept(FitInfo3{cross,fdout}.IndexMinMSE);
            prelabel=Test_dat'*B3{cross,fdout}(:,FitInfo3{cross,fdout}.IndexMinMSE) + FitInfo3{cross,fdout}.Intercept(FitInfo3{cross,fdout}.IndexMinMSE);
            A{cross,fdout}=B3{cross,fdout}(:,FitInfo3{cross,fdout}.IndexMinMSE);
            acrr(cross,fdout) = abs(corr( prelabel , Test_lab'));
%             train_SDI_cov = cov(Train_dat'*A{cross,fdout});
%             train_SDI_cov_norm = (A{cross,fdout}/train_SDI_cov)/(size(Train_dat,1)-1);
%             acti_patt(:,(cross-1)*10+fdout) = Train_dat*(Train_dat'*train_SDI_cov_norm );
            
            acti_patt(:,(cross-1)*10+fdout) = zscore(fast_haufe(Train_dat',A{cross,fdout},10));
%             cov(Train_dat')*A{cross,fdout}*inv(A{cross,fdout}'*Train_dat);
%             
%             (Train_dat'*train_SDI_cov_norm );
%             
            
%             train_cog = cov(Train_lab_new);
%             acti_patt(:,(cross-1)*10+fdout) = train_SDI_cov*A{cross,fdout}*inv(train_cog);
        end
    end
    
    
    precorr_A(out,:)=acrr(:);
    pre = precorr_A(out,:);
    pre_nan = pre(~isnan(pre));
    correlation(out) = mean(pre_nan)
    stdcorr(out) = std(pre_nan)
    acti_patt_nitem(:,:,out) = acti_patt;
end

% save('precorr_A','precorr_A')
% save('acti_patt_nitem','acti_patt_nitem')
% fast_haufe(train_dat',A{1,1})
