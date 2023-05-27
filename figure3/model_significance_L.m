
nIterNull=1%000;
load('data/behavioral_data'); 
load('data/SDI');
load('data/id_use.mat')
behav=behavioral_data(:,[13,14,27,28,30]);
array_null = zeros(size(behav,2),nIterNull);

parpool('local',8); %8核:local

tic
for i =1:size(behav,2)
    beha = behav(:,i);
    
    parfor j=1:nIterNull
        Posp = randperm(size(behavioral_data,1))';
        behap=beha(Posp);
        [train_set,test_set]= train_test(10,10);

        %% predict behavioral scores
        nSubj=877;nROI=360;
%         behav=behavioral_data(:,[13,14,27,28,30]);%cog
        
        nItem=size(behav,2);
        precorr=zeros(nItem,100);
        
        label = xlsread('data/glasser360_7networks.xlsx');
        label = label(:,3);
        
        
        vis = SDI(find(label==1),:);
        sm = SDI(find(label==2),:);
        da = SDI(label==3,:);
        va = SDI(find(label==4),:);
        lim = SDI(find(label==5),:);
        fp = SDI(find(label==6),:);
        dmn = SDI(find(label==7),:);
        
        
        % data=SDI;
        % data=[da;va;lim;fp;dmn];
         data=[vis;sm];
        
        lab=behap';
        acrr=zeros(10,10);
        % kfoldin=3;
        for cross=1%:10
            % kfoldout=10; % 10 CV for outer loop
            for fdout=1%:10
                fprintf('行为数据%d，置换检验%d，第%d次10折交叉验证%d折\n',i,j, cross,fdout);
                Train_dat = data(:,train_set{cross,fdout});
                Train_lab = lab(train_set{cross,fdout});
                Test_dat = data(:,test_set{cross,fdout});
                Test_lab = lab(test_set{cross,fdout});
                Train_lab_new=Train_lab;
                Test_lab_new=Test_lab;
                
                % % Test based on the optimal parameters obtained from the training data.
                
                %       indpout=1:size(Train_dat,1);
                %       [~,p_tmpout] = corr(Train_dat', Train_lab_new');
                %      indexout = find(p_tmpout<0.05);
                %      Train_dat = Train_dat(indexout,:);
                %      Test_dat = Test_dat(indexout,:);
                
                [B3,FitInfo3]=lasso(Train_dat',Train_lab_new','Alpha',0.01,'CV',3);
                prelabel=Test_dat'*B3(:,FitInfo3.IndexMinMSE) + FitInfo3.Intercept(FitInfo3.IndexMinMSE);
                A=B3(:,FitInfo3.IndexMinMSE);
                acrr(cross,fdout) = abs(corr( prelabel , Test_lab_new'));
            end
        end
        
        precorr(i,:)=acrr(:);
        pre = precorr(i,:);
        pre_nan = pre(~isnan(pre));
        correlation_perm_L(i,j) = mean(pre_nan);
        stdcorr(i,j) = std(pre_nan);
    end
%         array_null_L(i,j) = correlation;
end

toc

delete(gcp('nocreate'));