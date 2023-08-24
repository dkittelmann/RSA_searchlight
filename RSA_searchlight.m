%% RREPRESENTATIONAL SIMILARITY ANALYSIS 
% 
% Denise Kittelmann 


% This script performs a searchlight-based representational similarity analysis (RSA) 
% over the whole brain based on similarity template provided by TDT. 
% For each searchlight (i.e each voxel), a neural-based (here based on fMRI data) r
% epresentational dissimilarity matrix (RDM) is computed.
% Then, each neural RDM is correlated with two different stimulus-based RDMs.
% Finally, a brain mask (consisting of correlation coefficients in each voxel) 
% is created, which can be used later to extract the regions underlying 
% each searchlight.


% =====================================================================
%                           SET UP 
% =====================================================================


% REQUIREMENTS: 

% The Decoding Toolbox (TDT)*
% SPM12* 
% Betas from an SPM.mat (i.e. a GLM should have been estimated beforehand) 
% Model-based RDMs (e.g. stimuli-based RDMs)

% *both toolboxes need to be added to the MATLAB path!

clearvars; 


%% Create stimulus-based RDMs 

% within categories model: (high dissimilarity between stimulation & mental imagery, 
% low dissimilarity within stimulation & imagery conditions)
within_categories_model = [0 0 0 1 1 1]; 
within_categories_model_matrix = squareform(pdist(within_categories_model'));
%figure; imagesc(within_categories_model_matrix)


% stimulus type model: 
stimulus_type_model = [0 1 1 0.5 1 1; ...
                           1 0 1 1 0.5 1; ...
                           1 1 0 1 1 0.5; ...
                           0.5 1 1 0 1 1; ...
                           1 0.5 1 1 0 1; ...
                           1 1 0.5 1 1 0];
stimulus_type_model_matrix = stimulus_type_model;
%figure; imagesc(stimulus_type_model_matrix)


% Figure for all distance matrices 

x_labels = {'StimPress', 'StimFlutt', 'StimVibro', 'ImagPress', 'ImagFlutt', 'ImagVibro'}; 
y_labels = {'StimPress', 'StimFlutt', 'StimVibro', 'ImagPress', 'ImagFlutt', 'ImagVibro'};


RDMs = figure; 
subplot(1,2,1)
imagesc(within_categories_model_matrix)
title('Within Category RDM')
xticklabels(x_labels)
yticklabels(y_labels)
c= colorbar; 
c.Label.String = 'Dissimilarity'; 
subplot(1,2,2)
imagesc(stimulus_type_model_matrix)
title('Stimulus Type RDM') 
xticklabels(x_labels)
yticklabels(y_labels)
c= colorbar; 
c.Label.String = 'Dissimilarity'; 

% savesas(RDMs,fullfile(results_dir, 'RDMs.png'))

%% Compute fMRI-based RDMs for each sphere 

% searchlight radius is kept to default, i.e. 

for subj = 1:10 
    subj_str = sprintf('%02d', subj);

    % Set tdt decoding defaults
    cfg = decoding_defaults;
    cfg.analysis = 'searchlight';
    cfg.result.write = 1; 
    cfg.results.overwrite = 1; % set to 0 if you don't want to overwrite  

    % Set the results directory 
    cfg.results.dir = ['/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_', subj_str, '/'];

    % Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
    beta_loc = ['/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/data/sub_',  subj_str, '/1st_level_good_bad_Imag/'];

    % Set the filename of your brain mask (or your ROI masks as cell
    % matrix), e.g.'c:\exp\roi\roimaskright.img'
    cfg.files.mask = ['/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/data/sub_',  subj_str, '/1st_level_good_bad_Imag/mask.nii']; 

    % Set the label names to the regressor names for similarity analysis
    % e.g. 'button left' and 'button right'
    % Don't remember the names? -> run display_regressor_names(beta_loc)
    labelnames = {'StimPress', 'StimFlutt', 'StimVibro', 'ImagPress', 'ImagFlutt', 'ImagVibro'};

    % Set arbitrary labels 
    labels=[1 2 3 4 5 6];

    % set everything to dis-/similarity analysis (for available options as model parameters, check decoding_software/pattern_similarity/pattern_similarity.m)
    cfg.decoding.software = 'distance'; % calculates the distances instead of similarities
    cfg.decoding.method = 'classification';
    cfg.decoding.train.classification.model_parameters = 'cveuclidean2'; % cross validated euclidean
    cfg.results.output = 'other_average'; 
    % 'other_average' = average means across folds
    % 'other' = one similarity estimate per condition per run
    % 'other_meandist' = averages across (dis)similarity matrices of each cross-validation 
    %  iteration and across all cells of the lower diagonal (i.e. all distance comparisons)

    % Set additional parameters manually 

    %cfg.searchlight.unit = 'mm';
    %cfg.searchlight.radius = 12; 
    %cfg.searchlight.spherical = 1; if SL should be a spherical in mm (not in voxels)
    %cfg.verbose = 2; % you want all information to be printed on screen
    % cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q';
    % cfg.results.output = {'accuracy_minus_chance','AUC_minus_chance'};
    %cfg.results.output = {'accuracy_minus_chance'};

  
    % Nothing needs to be changed below for a standard similarity analysis using all data

    % The following function extracts all beta names and corresponding run
    % numbers from the SPM.mat
    regressor_names = design_from_spm(beta_loc);

    % Extract all information for the cfg.files structure (labels will be [1 2 3 4 5 6])
    cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_loc);

    % This creates a RSA design with cross-validation 
    cfg.design = make_design_similarity_cv(cfg);
   
    % Run decoding
    results = decoding(cfg);

end

% Plot RDM
% figure; heatmap(results.other_average.output{1}(1:6,1:6), 'colormap', jet)


%% Compute representational similiarity by correlating stimulus-based RDMs & fMRI based RDM 

 % Correlate sRDM with fRDM for each searchlight peak voxel for all
 % participants

 for subj = 1:10
     subj_str = sprintf('%02d', subj);

     results_dir = ['/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/sub_', subj_str, '/'];
     load(fullfile(results_dir,'res_other_average.mat'));
     mask_file = ['/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/data/sub_',  subj_str, '/1st_level_good_bad_Imag/mask.nii'];

     RDM_per_searchlight = results.other_average.output; % fRDM in each each sphere based on cv euclidean distance

     % Calculate correlation coefficients for each stimulus-based RDM

     RS_withincategoryRDMxnRDM = [];
     RS_stimulustypeRDMxnRDM = [];

     for i = 1:length(RDM_per_searchlight)
         RS_withincategoryRDMxnRDM(i) = corr(within_categories_model_matrix(:), RDM_per_searchlight{i}(:),'type','Spearman'); % spearman correlation
         RS_stimulustypeRDMxnRDM(i) = corr(stimulus_type_model_matrix(:), RDM_per_searchlight{i}(:),'type','Spearman');
     end

     % transpose each RSM to get a column vector
     RS_withincategoryRDMxnRDM = RS_withincategoryRDMxnRDM';
     RS_stimulustypeRDMxnRDM =  RS_stimulustypeRDMxnRDM';

     % Create brain mask with correlation coeffecients in each voxel for each participant

     % brain mask withincategory RDM x nRDM
     [Y,xyz] = spm_read_vols(spm_vol(mask_file));
     Y(results.mask_index)= RS_withincategoryRDMxnRDM;
     V = spm_vol(fullfile(results_dir,'beta_0001.nii')); % use header info from beta_001
     spm_write_vol(V,Y); % write mask

     % brain mask stimulus typ RDM x nRDM

     [Y,xyz] = spm_read_vols(spm_vol(mask_file));
     Y(results.mask_index)=  RS_stimulustypeRDMxnRDM;
     V = spm_vol(fullfile(results_dir,'beta_0002.nii')); % use header info from beta_002
     V = spm_write_vol(V,Y); % write mask

     save(fullfile(results_dir,'RS_withincategoryRDMxnRDM.mat'),'RS_withincategoryRDMxnRDM');
     save(fullfile(results_dir, 'RS_stimulustypeRDMxnRDM.mat'),'RS_stimulustypeRDMxnRDM');
 end


