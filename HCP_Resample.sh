#!/bin/bash
# PRE-PROCESSING DATA THEN RESAMPLE DATA ONTO HCP MESH AND FIND CONNECTIVITY
#-----------------------------------------------------------------------------#

ica=_ICA.ica
zero=0
wb_command=/home/zhallgri/workbench/bin_linux64/wb_command

SUBJECTS_DIR=$PWD

for SUBJECT in $1; do
	cd $SUBJECTS_DIR
	
	# Get run numbers for file hierarchy
	run_func=`cat ./$SUBJECT/scan_information.txt | grep RS_fMRI | awk '{print $1}'`

	echo $SUBJECT
	# Check if subject file na wb_command=/home/zhallgri/workbench/bin_linux64/wb_commandme folder is here -> for pre-processing the fmIR data to get hemispheres
	if [ ! -e ./$SUBJECT/subjectname ]; then
		touch ./$SUBJECT/subjectname
		echo $SUBJECT/struct_recon-all >> ./$SUBJECT/subjectname
	fi

	#Create directory to put de-noised processing in 
	mv ./$SUBJECT/bold ./$SUBJECT/unclean_fMRI
	mkdir ./$SUBJECT/bold ./$SUBJECT/bold/$zero$run_func

	#Put filtered functional in folder and unzip and re-name to match freesurfer hierarchy
	cp ./$SUBJECT/$SUBJECT$ica/filtered_func_data_clean.nii.gz ./$SUBJECT/bold/$zero$run_func/filtered_func_data_clean.nii.gz
	gunzip ./$SUBJECT/bold/*/filtered_func_data_clean.nii.gz 
	mv ./$SUBJECT/bold/*/filtered_func_data_clean.nii ./$SUBJECT/bold/$zero$run_func/f.nii

	# Pre-process the cleaned data to get hemispheres
	preproc-sess -s $SUBJECT -fsd bold -surface fsaverage lhrh -mni305 -fwhm 0 -nosmooth -per-run 	

	mkdir ./$SUBJECT/resampled-HCP
	
	# START RESAMPLING
	echo START RESAMPLING
	
	# Convert to gifti - WILL NEED TO MAKE THIS MORE GENERAL LATER
	gunzip ./$SUBJECT/bold/$zero$run_func/fmcpr.sm0.fsaverage.lh.nii.gz
	gunzip ./$SUBJECT/bold/$zero$run_func/fmcpr.sm0.fsaverage.rh.nii.gz 
	mri_convert ./$SUBJECT/bold/$zero$run_func/fmcpr.sm0.fsaverage.lh.nii ./$SUBJECT/resampled-HCP/lh_fsaverage.gii
	mri_convert ./$SUBJECT/bold/$zero$run_func/fmcpr.sm0.fsaverage.rh.nii ./$SUBJECT/resampled-HCP/rh_fsaverage.gii

	# NEED TO CHANGE /INGRIDLAB/ZHALLGRIMON WHEN MOVE INTO LOCAL FOLDERS
	# Resample left hemisphere
	$wb_command -metric-resample /ingridlab/zhallgrimson/RS/$SUBJECT/resampled-HCP/lh_fsaverage.gii \
		./standard_mesh_atlases/resample_fsaverage/fsaverage_std_sphere.L.164k_fsavg_L.surf.gii \
		./standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii \
		ADAP_BARY_AREA \
		./$SUBJECT/resampled-HCP/resampled.surface.lh.32k_fs_LR.func.gii \
		-area-metrics \
		./standard_mesh_atlases/resample_fsaverage/fsaverage.L.midthickness_va_avg.164k_fsavg_L.shape.gii \
		./standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii 

	# Resample right hemisphere
	 $wb_command -metric-resample /ingridlab/zhallgrimson/RS/$SUBJECT/resampled-HCP/rh_fsaverage.gii \
		./standard_mesh_atlases/resample_fsaverage/fsaverage_std_sphere.R.164k_fsavg_R.surf.gii \
		./standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii \
		ADAP_BARY_AREA \
		./$SUBJECT/resampled-HCP/resampled.surface.rh.32k_fs_LR.func.gii \
		-area-metrics \
		./standard_mesh_atlases/resample_fsaverage/fsaverage.R.midthickness_va_avg.164k_fsavg_R.shape.gii \
		./standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii
	
	cd ./$SUBJECT/resampled-HCP
	
	# Create Dense Time series
	$wb_command -cifti-create-dense-timeseries resampled_series.dtseries.nii -left-metric resampled.surface.lh.32k_fs_LR.func.gii -right-metric resampled.surface.rh.32k_fs_LR.func.gii

	# START CONNECTIVITY ANALYSIS
	echo START CONNECTIVITY

	# Create standard normal distribution 
	$wb_command -cifti-reduce resampled_series.dtseries.nii MEAN mean.dscalar.nii   
	$wb_command -cifti-reduce resampled_series.dtseries.nii STDEV stdev.dscalar.nii
	$wb_command -cifti-math "(x - mean) / stdev" norm_resampled_series.dtseries.nii -fixnan 0 -var x resampled_series.dtseries.nii -var mean mean.dscalar.nii -select 1 1 -repeat -var stdev stdev.dscalar.nii -select 1 1 -repeat
 
	# Parcellate dtseries
	$wb_command -cifti-parcellate norm_resampled_series.dtseries.nii /ingridlab/zhallgrimson/RS/parcellation-180.dlabel.nii COLUMN parcel_norm_resampled_series.ptseries.nii -method MEAN

	# Correlate parcellation file 
	end=_connect.pconn.nii
	$wb_command -cifti-correlation  parcel_norm_resampled_series.ptseries.nii $SUBJECT$end

done

echo Finished 
