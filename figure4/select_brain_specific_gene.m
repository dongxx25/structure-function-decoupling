clear all;
clc;

gene_all  = xlsread('data/specific_gene.xlsx');
gene_spe = readtable('data/specific_gene.xlsx');
gene_id = table2array(gene_spe(:,1));

% gene_cate = table2array(gene_spe(:,29));
% idx = find(ismember(gene_cate, 'Not detected' ));
% gene_all(idx,:) = [];
% gene_id(idx,:) = [];

gene_brain = gene_all(:,7);

% gene_all(:,[7])=[];

gene_all_med = 5 * median(transpose(gene_all));

% sum(gene_brain > gene_all_med)


gene_brain_speci = [];
pos_gene_brain_speci = [];

for i = 1:size(gene_all,1)
    
    if gene_brain(i) > gene_all_med(i)
        gene_brain_new = gene_brain(i);
        pos = i;
        gene_brain_speci = [gene_brain_speci;gene_brain_new];
        pos_gene_brain_speci = [pos_gene_brain_speci;pos];
    end
    
end

gene_brain_specific = gene_id(pos_gene_brain_speci);

%  save('gene_brain_specific','gene_brain_specific')





