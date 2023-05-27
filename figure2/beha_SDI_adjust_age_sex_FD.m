
clear all;
clc;

%% Load some tools
 addpath(genpath('PALM')); %https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM
 addpath(genpath('FSLNets')); %https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets

 load('/data7/dxx/SDI/full_code/code/data/beha.mat')
 load('/data7/dxx/SDI/full_code/code/data/SDI.mat')
 behavioral_data = beha;
 
 SDI = array2table(SDI');
 SDI = [SDI,behavioral_data(:,[111,112,113])];
%% Adjust behavioral data for age and sex
behavioral_data_resid = nan((height(behavioral_data)),(width(behavioral_data)-4));
for i = 1:size(behavioral_data_resid,2)    
    Yname = cell2mat(behavioral_data.Properties.VariableNames(i));
    modelspec = [Yname ' ~ Gender + Age_in_Yrs '];
    lm = fitlm(behavioral_data,modelspec);
    behavioral_data_resid(:,i)=lm.Residuals.Raw;
end
% save('new/behavioral_data_resid','behavioral_data_resid')

% behavioral_data_resid=behavioral_data_resid(:,[4:size(behavioral_data_resid,2)]);

%% Adjust SDI for FD
SDI_residFD = nan(height(SDI),width(SDI)-3);

for i = 1:size(SDI_residFD,2)    
    Yname = cell2mat(SDI.Properties.VariableNames(i));
    modelspec = [Yname ' ~ FD '];
    lm = fitlm(SDI,modelspec);
    SDI_residFD(:,i)=lm.Residuals.Raw;
end

%% Adjust SDI for age sex & FD
SDI_resid_age_sex_FD = nan(height(SDI),width(SDI)-3);

for i = 1:size(SDI_resid_age_sex_FD,2)    
    Yname = cell2mat(SDI.Properties.VariableNames(i));
    modelspec = [Yname ' ~Gender + Age_in_Yrs + FD '];
    lm = fitlm(SDI,modelspec);
    SDI_resid_age_sex_FD(:,i)=lm.Residuals.Raw;
end

%save('new/SDI_resid','SDI_resid')


%% without adjust for behaviral data
behavioral_data = table2array(behavioral_data);
behavioral_data = behavioral_data(:,(1:109));




%% normalize
SDI_residFD = nets_demean(SDI_residFD);
SDI_residFD = SDI_residFD./nanstd(SDI_residFD);

SDI_resid_age_sex_FD = nets_demean(SDI_resid_age_sex_FD);
SDI_resid_age_sex_FD = SDI_resid_age_sex_FD./nanstd(SDI_resid_age_sex_FD);

behavioral_data_resid = nets_demean(behavioral_data_resid);
behavioral_data_resid = behavioral_data_resid./nanstd(behavioral_data_resid);

behavioral_data_resid = palm_inormal(behavioral_data_resid);

behavioral_data = nets_demean(behavioral_data);
behavioral_data = behavioral_data./nanstd(behavioral_data);

behavioral_data = palm_inormal(behavioral_data);


save('new/behavioral_data_resid','behavioral_data_resid')
save('new/behavioral_data','behavioral_data')
save('new/SDI_residFD','SDI_residFD')
save('new/SDI_resid_age_sex_FD','SDI_resid_age_sex_FD')
