% Average decoding mean correlation 
%
% Denise Kittelmann 

% This script calculates the mean correlation/decoding accuracy for each
% cluster

% ------------------------------------------------------------------------

resultsfig_dir = '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/figures/';


%% Calculate average mean accruacy for  across category model in each significant cluster

mean_cluster_acc1 = []; 

for i= 1:5

% read in averaged correlation maps across participants (calculated im imCalc in spm12) & the significant cluster map from the 2nd level analysis
averagedcormap = '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/res_betweenconditions/corabetween_averaged_across_participants.nii';
clustermask = ['/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/res_betweenconditions/cluster', num2str(i), 'binary.nii']; 
% binary mask for each significant clusters were extracted from spm results


hdr = spm_vol(averagedcormap); % get hdr info from averaged correlation map
vol = spm_read_vols(hdr); % read into volumne
maskvol = spm_read_vols(spm_vol(clustermask));
cluster_acc = vol(logical(maskvol));
mean_cluster_acc1(i) = mean(cluster_acc); % calculate mean correlation 

end 


%% Calculate average mean accruacy for the within category model in each significant cluster

mean_cluster_acc2 = []; 

for i= 1:7
averagedcormap = '/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/res_stimulustype/coracross_averaged_within_participants.nii';
clustermask = ['/Users/denisekittelmann/Documents/MATLAB/RSA_statsproject/results/res_stimulustype/cluster', num2str(i), 'binary.nii'];


hdr = spm_vol(averagedcormap); % get hdr info from averaged correlation map
vol = spm_read_vols(hdr); % read into volumne
maskvol = spm_read_vols(spm_vol(clustermask));
cluster_acc = vol(logical(maskvol));
mean_cluster_acc2(i)= mean(cluster_acc); % calculate mean correlation 
end

%% Plot average correlations


cluster_labels = {'rSII', 'lSII', 'rSI', 'lIFG', 'lAin'}; 

meancor = figure; bar(mean_cluster_acc1)
xticklabels(cluster_labels); 
ylabel('mean correlation'); 
ylim([0,1])


saveas(meancor,fullfile(resultsfig_dir, 'meancorrelation.png'))



