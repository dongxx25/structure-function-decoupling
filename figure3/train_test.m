        
function [train_set,test_set]= train_test(i,j)

% split train set and test set
        load('data/id_use.mat')
        ID = readtable('data/family_ID.xlsx');
        ID_all = table2array(ID(:,1));
        [c, ia, ib] = intersect(ID_all, id_use);
        
        Family_ID = table2cell(ID(:,2));
        Family_ID = Family_ID(ia);
        nSubj=size(Family_ID,1);
        train_set=cell(i,j);
        test_set=cell(i,j);
        Famtrain_set=cell(i,j);
        Famtest_set=cell(i,j);
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
        