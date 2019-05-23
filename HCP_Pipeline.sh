#!/bin/bash
# PRE-PROCESS HCP SUBJECT DATA THROUGH FSL AND FREESURFER
# ----------------------------------------------------------
# Input: must be in form melodic_pp "SUBJECT1 SUBJECT2 SUBJECT3 SUBJECT4" <FIX>
#	-"SUBJECT1 SUBJECT2 SUBJECT3 SUBJECT4" is list of subjects that can be processed
#	-<FIX>: if want to use fix to extract components then put in: classify  
#	-<FIX>: if want to classify and removal then put in: remove

# If nothing passed assume skip fix
var=${2:skip}

# Directory with subjects
SUBJECTS_DIR=/ingridlab/mark/HCP_Subject 
cd /ingridlab/mark/HCP_Subject

# Used for batch 
end=_ICA
mel=.ica
tmp=_struct_tmp
act=struct_recon-all
boldpathend='"'

# Loop through group of subjects
for SUBJECT in $1; do 
	echo $SUBJECT	

	# SKULL STRIP
	bet ./$SUBJECT/Struct/MNINonLinear/T1w_restore.1.60.nii.gz ./$SUBJECT/Struct/MNINonLinear/T1w_restore.1.60_brain.nii.gz

	# Edit design file with subject name 
	cp preprocess_design_HCP.fsf ./$SUBJECT/preprocess_design_${SUBJECT}.fsf
	# Find path to function
	sed -i -e "s/SUBJECT/$SUBJECT/" ./$SUBJECT/preprocess_design_${SUBJECT}.fsf #replace names
	sed -i -e "s/SUBJECT_7T_rfMRI_REST1_PA.nii.gz/${SUBJECT}_7T_rfMRI_REST1_PA.nii.gz/" ./$SUBJECT/preprocess_design_${SUBJECT}.fsf #replace names
	sed -i -e "s/SUBJECT_ICA/${SUBJECT}${end}/" ./$SUBJECT/preprocess_design_${SUBJECT}.fsf #replace ica path
	
	# MELODIC
	# motion correction, co-registration, high pass filtering, and ICA
	feat ./$SUBJECT/preprocess_design_${SUBJECT}.fsf 
	
	# FIX 
        #if [ "$2" = classify ]; then 
	#	~/Documents/fix/fix -f ./$SUBJECT/$SUBJECT$end$mel
	# Change TRAIN.RData
	#	~/Documents/fix/fix -c ./$SUBJECT/$SUBJECT$end$mel TRAIN.RData 20
	#elif [ "$2" = remove ]; then
	#	~/Documents/fix/fix ./$SUBJECT/$SUBJECT$end$mel TRAIN.RData 20
	#fi

done 

