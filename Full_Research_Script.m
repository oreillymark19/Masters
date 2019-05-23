%% Main Analysis Script

%% Read in data
% Add necessary packages to path
addpath('/mnt/Data/gifti-1.6/');
glasser_index_check = dlmread('reg_index.csv');

% Get specifics from user on what atlas should be used
prompt = 'What parcellation would you like to use, glasser, yeo, or ICA? ';
atlas{1} = input(prompt,'s');
atlas{1} = lower(atlas{1});
atlas{1} = atlas{1}(~isspace(atlas{1}));
fprintf('\n');

switch atlas{1}
    % Load in proper atlas based on user input
    case 'yeo'
        parcel_prompt = 'How many parcels will be used, 100 200 300 400 500 600 800 1000? ';
        atlas{2} = input(parcel_prompt);
        fprintf('\n')
        network_num_prompt =  'How many networks, 7 or 17? ';
        atlas{3} = input(network_num_prompt);
        fprintf('\n')
        addpath('/mnt/Data/fieldtrip-20180312');
        yeo = ft_read_cifti(['/mnt/Data/git-repos/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_', ...
            num2str(atlas{2}) ,'Parcels_',num2str(atlas{3}) ,'Networks_order.dlabel.nii'],'mapname','array');
        
        labels = yeo.dlabel;
        if atlas{3} == 7
                    yeo_networks = ["Vis", 
                        "SomMot", 
                        "DorsAttn",
                        "SalVentAttn",
                        "Limbic",
                        "Cont",
                        "Default"];
                elseif atlas{3} == 17
                    yeo_networks = ["VisCent",
                        "VisPeri",
                        "SomMotA",
                        "SomMotB",
                        "DorsAttnA", 
                        "DorsAttnB",
                        "SalVentAttnA",
                        "SalVentAttnB",
                        "Limbic_OFC",
                        "Limbic_TempPole",
                        "ContA",
                        "ContB",
                        "ContC",
                        "DefaultA",
                        "DefaultB",
                        "DefaultC",
                        "TempPar"]';
                end
                
                for i = 1:size(yeo_networks',2)
                    yeo_network = yeo_networks(i);
                    temp_mask = contains(yeo.dlabellabel, yeo_networks(i));
                    mask = find(temp_mask == 1);
                    for j = 1:size(mask,2)
                        index = mask(j);
                        labels(find(labels(:,1) == index),2) = i;
                    end
                    clear i temp_mask mask yeo_network index
                end
                
                load('bad_index_list.mat')
                labels(bad_indices,:) = [];
        
    case 'glasser'
        glasser = gifti(['/mnt/Data/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Glasser_Parcellation32k_fs_LR.dlabel.gii']);
        labels = glasser.cdata;
        index_check = dlmread('reg_index.csv');
        
    case 'ica'
        ICA_exist_prompt = 'Has MELODIC already been run on the data, Y or N? ';
        ICA_exist = input(ICA_exist_prompt, 's');
        ICA_exist = lower(ICA_exist);
        ICA_exist = ICA_exist(~isspace(ICA_exist));
        fprintf('\n')
    
        switch ICA_exist
            case 'y'
                mov_or_rest_prompt = 'Would you like to use the rest [R] or movie [M] components? '
                mov_or_rest = input(mov_or_rest_prompt, 's');
                mov_or_rest = lower(mov_or_rest);
                mov_or_rest = mov_or_rest(~isspace(mov_or_rest));
                fprintf('\n')
                
                switch mov_or_rest
                    case 'm'
                        ICA = gifti(['/mnt/Data/New_Subjects/Full_4mm_smoothed_ICA_Movie/melodic_IC.gii']);
                        ICA_SM = ICA.cdata';
                        
                    case 'r'
                        ICA = gifti(['/mnt/Data/New_Subjects/Full_4mm_smoothed_ICA_Rest/melodic_IC.gii']);
                        ICA_SM = ICA.cdata';
                        
                    otherwise
                        fprintf('Error: You entered an invalid response')
                end
                
            case 'n'
                fprintf('You need to run MELODIC first \n\n')
                return
                
            otherwise
                fprintf('Error: You entered an invalid response')
        end
        
    otherwise
        fprintf('Error: You did not enter a valid atlas. \n\n')
        
end

clear parcel_prompt network_num_prompt prompt

% Ask is user wants to read in data
prompt = 'Would you like to read in data, Y or N? [N]: ';
read_in = input(prompt,'s');
read_in = lower(read_in);
if isempty(read_in)
    read_in = 'n'
end
clear prompt
fprintf('\n')

switch read_in
    
    case 'n'
        fprintf('No data was read \n\n');
        
    case 'y'
        list_prompt = 'Name of subjects list to use: ';
        sub_list = input(list_prompt,'s');
        filename = ['/mnt/Data/Subject_Lists/', sub_list];
        fidi = fopen(filename, 'rt');
        D = textscan(fidi, '%s%s', 'CollectOutput',1);
        subjects_list = D{1,1};
        fprintf('\n')
        clear D fidi filename list_prompt
        
        run_prompt = 'Which runs would you like to read in? If more than one, enter in square brackets, i.e. [1 2 3] ';
        runs = input(run_prompt);
        if runs > 4;
            disp('That is too many, max is 4')
            return
        end
        
        % Declare variable to hold subjects with missing data
        missing_data = [];
        count = 0;
        
        switch atlas{1}
            case {'yeo', 'glasser'}
                
                % Loop through subject list
                for i = 1:size(subjects_list,1)
                    subject = subjects_list{i,1};
                    fprintf([subject, '\n']);
                    filepath = ['/mnt/Data/New_Subjects/' ,num2str(subject),'/MNINonLinear/Results/'];
                    for r = 1:size(runs,2)
                        run = runs(r);
                        if run == 1
                            rest_dir = 'PA';
                            mov_dir = 'AP';
                        elseif run == 2
                            rest_dir = 'AP';
                            mov_dir = 'PA';
                        elseif run == 3
                            rest_dir = 'PA';
                            mov_dir = 'PA';
                        elseif run == 4
                            rest_dir = 'AP';
                            mov_dir = 'AP';
                        end
                        
                        % Check if subject files exist
                        if isfile([filepath,'rfMRI_REST', num2str(run),'_7T_',rest_dir,'/rfMRI_REST',num2str(run),'_7T_',rest_dir,'_Atlas_MSMAll_hp2000_clean.dtseries.gii'])==0|...
                                isfile([filepath,'tfMRI_MOVIE', num2str(run),'_7T_',mov_dir,'/tfMRI_MOVIE', num2str(run),'_7T_',mov_dir,'_Atlas_MSMAll_hp2000_clean.dtseries.gii'])==0
                            fprintf('Subject %s is missing run %d \n', subject, run);
                            count = count + 1;
                            missing_data(count,1) = i;
                            missing_data(count,2) = run;
                            continue
                        end
                        
                        % Read in rest data for subject and remove subcortical data
                        rest_data = gifti([filepath, 'rfMRI_REST', num2str(run),'_7T_',rest_dir,'/rfMRI_REST',num2str(run),'_7T_',rest_dir,'_Atlas_MSMAll_hp2000_clean.dtseries.gii']);
                        rest_time_series = rest_data.cdata(1:64984,:);
                        
                        %Read in movie data for subject and remove subcortical data
                        movie_data = gifti([filepath,'tfMRI_MOVIE', num2str(run),'_7T_',mov_dir,'/tfMRI_MOVIE', num2str(run),'_7T_',mov_dir,'_Atlas_MSMAll_hp2000_clean.dtseries.gii']);
                        mov_time_series = movie_data.cdata(1:64984,:);
                        
                        % Loop through regions
                        for j = 1:max(labels)
                            % Isolate vertices corresponding to region j
                            mov_roi = mov_time_series(labels(:,1)==j,:);
                            rest_roi = rest_time_series(labels(:,1)==j,:);
                            
                            % Average vertices in region j to get single time series
                            mov_mean_roi(j,:,i,r) = mean(mov_roi(:,1:900));
                            rest_mean_roi(j,:,i,r) = mean(rest_roi(:,1:900));
                            
                            % Clearing useless variables
                            clear mov_roi rest_roi
                        end

                        % Clearing useless variables
                        clear movie_data rest_data ext mov_time_series rest_time_series
                    end
                    
                end
                clear i j subject count
                
                missing_subjects = unique(missing_data(:,1));
                
                mov_mean_roi(:,:,missing_subjects,:) = [];
                rest_mean_roi(:,:,missing_subjects,:) = [];
                
            case 'ica'
                % Read in subject data, apply centering, and apply dual regression
                fprintf('Applying Dual Regression \n\n')
                
                for i = 1:size(subjects_list)
                    % Load in subject number
                    subject = subjects_list{i,1};
                    fprintf([subject, '\n']);
                    filepath = ['/mnt/Data/New_Subjects/' ,num2str(subject),'/MNINonLinear/Results/'];
                    
                    % Check if subject files exist
                    if isfile([filepath,'rfMRI_REST2_7T_AP/rfMRI_REST2_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries.gii'])==0|...
                            isfile([filepath,'tfMRI_MOVIE2_7T_PA/tfMRI_MOVIE2_7T_PA_Atlas_MSMAll_hp2000_clean.dtseries.gii'])==0
                        fprintf('Subject %s is missing a file \n', subject);
                        count = count + 1;
                        missing_data(count) = i;
                        %continue
                    end
                    
                    % Load in subject data - Rest
                    rest_data = gifti([filepath, 'rfMRI_REST1_7T_PA/4mm_smoothed_rfMRI_REST1_MSMAll.dtseries.gii']);
                    rest_time_series = rest_data.cdata;
                    rest_time_series(any(isnan(rest_time_series), 2), :) = [];
                    
                    % Load in subject data - Movie
                    movie_data = gifti([filepath, 'tfMRI_MOVIE1_7T_AP/4mm_smoothed_tfMRI_MOVIE1_MSMAll.dtseries.gii']);
                    movie_time_series = movie_data.cdata(:, 1:size(rest_time_series,2));
                    movie_time_series(any(isnan(movie_time_series), 2), :) = [];
                    
                    % Centering
                    [voxels time] = size(rest_time_series);
                    components = size(ICA_SM,1);
                    
                    % Spatial Maps
                    cent_maps1 = (eye(voxels/2) - (1/(voxels))*ones(voxels/2,1)*ones(voxels/2,1)')*(ICA_SM(:,1:45641)');
                    cent_maps2 = (eye(voxels/2) - (1/(voxels))*ones(voxels/2,1)*ones(voxels/2,1)')*(ICA_SM(:,45642:91282)');
                    cent_maps = [cent_maps1; cent_maps2];
                    
                    clear cent_maps1 cent_maps2
                    
                    % Data - Time Centering
                    Y_time_rest = (eye(time) - (1/(time))*ones(time,1)*ones(time,1)')*(rest_time_series');
                    Y_time_movie = (eye(time) - (1/(time))*ones(time,1)*ones(time,1)')*(movie_time_series');
                    
                    % Data - Spatial Centering
                    Y_cent1 = (eye(voxels/2) - (1/(voxels))*ones(voxels/2,1)*ones(voxels/2,1)')*(Y_time_rest(:,1:45641)');
                    Y_cent2 = (eye(voxels/2) - (1/(voxels))*ones(voxels/2,1)*ones(voxels/2,1)')*(Y_time_rest(:,45642:91282)');
                    Y_cent_rest = [Y_cent1; Y_cent2];
                    
                    clear Y_cent1 Y_cent2 Y_time_rest
                    
                    Y_cent1 = (eye(voxels/2) - (1/(voxels))*ones(voxels/2,1)*ones(voxels/2,1)')*(Y_time_movie(:,1:45641)');
                    Y_cent2 = (eye(voxels/2) - (1/(voxels))*ones(voxels/2,1)*ones(voxels/2,1)')*(Y_time_movie(:,45642:91282)');
                    Y_cent_movie = [Y_cent1; Y_cent2];
                    
                    clear Y_cent1 Y_cent2 Y_time_movie
                    
                    % First Regression
                    rest_mean_roi(:,:,i) = linsolve(cent_maps,Y_cent_rest);
                    mov_mean_roi(:,:,i) = linsolve(cent_maps,Y_cent_movie);
                    
                    % Second Regression
                    S_rest_ind(:,:,i) = linsolve(rest_mean_roi(:,:,i)', Y_cent_rest');
                    S_movie_ind(:,:,i) = linsolve(rest_mean_roi(:,:,i)', Y_cent_movie');
                    
                    
                    clear Y_cent_movie Y_cent_rest cent_maps voxels time components rest_data rest_time_series i                   
                    
                end
                
                clear i j subject count
                mov_mean_roi(:,:,missing_data) = [];
                rest_mean_roi(:,:,missing_data) = [];
                
                rest_mean_roi = permute(rest_mean_roi, [2,1,3]);
                mov_mean_roi = permute(mov_mean_roi, [2 1 3]);
                
        end
        
    otherwise
        fprintf('Error: Your entry is invalid \n\n')
        
end

fprintf('\n')
clear read_in

%% Analysis Techniques

%% Variance
prompt = 'Would you like to use Fisher transformed values, Y or N? [N]: ';
fish = input(prompt,'s');
fish = lower(fish);
fish = fish(~isspace(fish));
clear prompt
fprintf('\n')
        
if strcmp(fish, 'y')
    z_vals = 1;
else
    z_vals = 0;
end
clear fish

analysis_prompt = 'Which analysis technique would you like to use? (Full Matrix, Within Network, Fingerprint): ';
analysis_choice = input(analysis_prompt, 's');
analysis_choice = lower(analysis_choice);
analysis_choice = analysis_choice(~isspace(analysis_choice));
fprintf('\n\n')

switch  analysis_choice
    case 'fullmatrix'
        fprintf('Full Matrix \n\n')
                [rest_var, rest_corr, rest_range_var, rest_range_corr] = full_mat_var(rest_mean_roi, 10, z_vals);
                [movie_var, movie_corr, movie_range_var, movie_range_corr] = full_mat_var(mov_mean_roi, 10, z_vals);                             
        
    case 'withinnetwork'
        fprintf('Within Network \n\n')
        prompt = 'Which network(s) would you like to examine? If more than one, enter in square brackets, i.e. [1 2 3] ';
        networks = input(prompt);
        clear prompt
        fprintf('\n')
        
        if any(networks > 7|networks < 1) && strcmp(atlas{1}, 'yeo')
            fprintf('Error: You entered a network number outside the bounds of possible values \n\n')
            return
        elseif any(networks > max(glasser_index_check(:,2))|networks < 1) && strcmp(atlas{1}, 'glasser')
            fprintf('Error: You entered a network number outside the bounds of possible values \n\n')
            return    
        elseif any(networks > size(mov_mean_roi,1)|networks < 1) && strcmp(atlas{1}, 'ica')
            fprintf('Error: You entered a network number outside the bounds of possible values \n\n')
            return  
        end
        
        prompt = 'Would you like to examine these networks individually [Y] or together as a larger network [N]? ICA must be [N]';
        read_in = input(prompt,'s');
        read_in = lower(read_in);
        read_in = read_in(~isspace(read_in));
        if isempty(read_in)
            read_in = 'y';
        end
        clear prompt
        fprintf('\n')
        
        switch read_in
            case 'y'
                for i = 1:size(networks,2)
                    switch atlas{1}
                        case 'yeo'
                            network = unique(labels(labels(:,2) == networks(i),1));
                        case 'glasser'
                            network = glasser_index_check(glasser_index_check(:,2) == networks(i),1);
                        case 'ica'
                            fprintf('Error: Not possible to analyze single networks in ICA \n\n')
                            return
                    end
                    [wn_rest_var(:,i) wn_rest_corr(:,i) wn_rest_range_var(:,i), wn_rest_range_corr(:,i)] = within_network_var(rest_mean_roi, 10, network, z_vals);
                    [wn_movie_var(:,i) wn_movie_corr(:,i) wn_movie_range_var(:,i), wn_movie_range_corr(:,i)] = within_network_var(mov_mean_roi, 10, network, z_vals);
                    clear network
                end
                
            case 'n'
                switch atlas{1}
                    case 'yeo'
                        network = unique(labels(find(sum(labels(:,2) == networks,2)),1));
                    case 'glasser'
                        network = glasser_index_check(find(sum(glasser_index_check(:,2) == networks,2)==1),1);
                    case 'ica'
                        network = networks;
                end
                
                [wn_rest_var wn_rest_corr wn_rest_range_var, wn_rest_range_corr] = within_network_var(rest_mean_roi, 10, network, z_vals);
                [wn_movie_var wn_movie_corr wn_movie_range_var, wn_movie_range_corr] = within_network_var(mov_mean_roi, 10, network, z_vals);
                
            otherwise
                fprintf('Error: Your entry is invalid \n\n')
        end
        
    case 'fingerprint'
        fprintf('Fingerprint \n\n')
        prompt = 'Which network(s) would you like to examine? If more than one, enter in square brackets, i.e. [1 2 3] ';
        networks = input(prompt);
        clear prompt
        fprintf('\n')
        
        if any(networks > 7|networks < 1) && strcmp(atlas{1}, 'yeo')
            fprintf('Error: You entered a network number outside the bounds of possible values \n\n')
            return
        elseif any(networks > max(glasser_index_check(:,2))|networks < 1) && strcmp(atlas{1}, 'glasser')
            fprintf('Error: You entered a network number outside the bounds of possible values \n\n')
            return
        elseif any(networks > size(mov_mean_roi,1)|networks < 1) && strcmp(atlas{1}, 'ica')
            fprintf('Error: You entered a network number outside the bounds of possible values \n\n')
            return
        end
        
        prompt = 'Would you like to examine these networks individually [Y] or together as a larger network [N]? ';
        read_in = input(prompt,'s');
        read_in = lower(read_in);
        read_in = read_in(~isspace(read_in));
        if isempty(read_in)
            read_in = 'y';
        end
        clear prompt
        fprintf('\n')
        
        switch read_in
            case 'y'
                for i = 1:size(networks,2)
                    switch atlas{1}
                        case 'yeo'
                            network = unique(labels(labels(:,2) == networks(i),1));
                        case 'glasser'
                            network = glasser_index_check(glasser_index_check(:,2) == networks(i),1);
                        case 'ica'
                            network = networks(i);
                    end
                    [fp_rest_var(:,i) fp_rest_corr(:,i) fp_rest_range_var(:,i), fp_rest_range_corr(:,i)] = fingerprint_var(rest_mean_roi, 10, network, z_vals);
                    [fp_movie_var(:,i) fp_movie_corr(:,i) fp_movie_range_var(:,i), fp_movie_range_corr(:,i)] = fingerprint_var(mov_mean_roi, 10, network, z_vals);
                    clear network
                end
                
            case 'n'
                switch atlas{1}
                    case 'yeo'
                        network = unique(labels(find(sum(labels(:,2) == networks,2)),1));
                    case 'glasser'
                        network = glasser_index_check(find(sum(glasser_index_check(:,2) == networks,2)==1),1);
                    case 'ica'
                        network = networks;
                end
                
                [fp_rest_var fp_rest_corr fp_rest_range_var, fp_rest_range_corr] = fingerprint_var(rest_mean_roi, 10, network, z_vals);
                [fp_movie_var fp_movie_corr fp_movie_range_var, fp_movie_range_corr] = fingerprint_var(mov_mean_roi, 10, network, z_vals);
                
            otherwise
                fprintf('Error: Your entry is invalid \n\n')
        end
        
    otherwise
        fprintf('Error: Your entry is invalid \n\n')
            
end

clear read_in i

%% Clear Analysis

clearvars -except atlas filepath labels missing_data mov_mean_roi rest_mean_roi sub_list subjects_list yeo yeo_networks glasser glasser_index_check S_rest_ind S_movie_ind
clc


%% Yeo Individual Variability Paper
diff_null_dist = zeros(5000,1);
mov_null_dist = zeros(5000,1);
rest_null_dist = zeros(5000,1);

for p = 1:100
    for run = 1:size(mov_mean_roi,4)
        for i = 1:size(mov_mean_roi,3)
            surr_data_mov(:,:,i,run) = phaseran(mov_mean_roi(:,:,i,run)',1);
            surr_data_rest(:,:,i,run) = phaseran(rest_mean_roi(:,:,i,run)',1);
        end
    end
    
    [mov_var_region, mov_var_sess] = anti_corr(mov_mean_roi);
    [rest_var_region, rest_var_sess] = anti_corr(rest_mean_roi);

% Uncomment this to include tSNR in regression model
%     for run = 1:4
%         for roi = 1:360
%             for sub = 1:size(mov_mean_roi,3)
%                 mov_mean = mean(mov_mean_roi(roi,:,sub,run));
%                 rest_mean = mean(rest_mean_roi(roi,:,sub,run));
%                 mov_std = std(mov_mean_roi(roi,:,sub,run));
%                 rest_std = std(rest_mean_roi(roi,:,sub,run));
%                 
%                 mov_sub_tSNR(sub) = mov_mean/mov_std;
%                 rest_sub_tSNR(sub) = rest_mean/rest_std;
%                 
%                 clear mov_mean rest_mean mov_std rest_std
%             end
%             mov_roi_tSNR(roi,run) = mean(mov_sub_tSNR);
%             rest_roi_tSNR(roi,run)= mean(rest_sub_tSNR);
%             
%             clear mov_sub_tSNR rest_sub_tSNR
%         end
%     end
%     
%     mov_tSNR = mean(mov_roi_tSNR,2);
%     rest_tSNR = mean(rest_roi_tSNR,2);
%     
%     clear mov_roi_tSNR rest_roi_tSNR
%     
%     rest_var_sess(:,2) = rest_tSNR;
%     mov_var_ses(:,2) = mov_tSNR;
    
    %Add intercept term
    rest_var_sess(:,2) = ones(size(rest_var_sess,1),1);
    mov_var_sess(:,2) = ones(size(mov_var_sess,1),1);
    
    [rest_beta,rest_sigma,rest_resid] = mvregress(rest_var_sess,rest_var_region);
    [mov_beta,mov_sigma,mov_resid] = mvregress(mov_var_sess,mov_var_region);
    
    mean_mov_resid = mean(mov_resid,2);
    mean_rest_resid = mean(rest_resid,2);
    
    mov_null_dist(p) = median(mean_mov_resid);
    rest_null_dist(p) = median(mean_rest_resid);
    
    diff_reg_var = mean_rest_resid - mean_mov_resid;
    diff_null_dist(p) = median(diff_reg_var);
    
    clear diff_reg_var mean_mov_reg_var mean_rest_reg_var mov_reg_var rest_reg_var rest_beta rest_sigma rest_resid ...
        mov_beta mov_sigma mov_resid rest_var_sess mov_var_sess mov_var_region mov_var_sess rest_var_region rest_var_sess subs
    
end

mean_mov_var = mean(mov_var_region,2);
mean_rest_var = mean(rest_var_region,2);

diff_mean_var = mean_rest_var - mean_mov_var;

diff_reg_var = rest_reg_var - mov_reg_var;

% Without regression of intrasub variability

mean_mov_var_region = mean(mov_var_region,2);
mean_rest_var_region = mean(rest_var_region,2);

diff_var_region = mean_rest_var_region - mean_mov_var_region;

%% Binned DFC Analysis

% Full Matrix Analysis
z_vals = 0;
mov_var_bins = binned_DFC(mov_mean_roi, z_vals)
rest_var_bins = binned_DFC(rest_mean_roi, z_vals)

mean_mov_bins = mean(mov_var_bins,3);
mean_rest_bins = mean(rest_var_bins,3);





   
