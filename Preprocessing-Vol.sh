#!/bin/bash
# PRE-PROCESS SUBJECT DATA THROUGH FSL AND FREESURFER
# ----------------------------------------------------------
# Input: must be in form melodic_pp "SUBJECT1 SUBJECT2 SUBJECT3 SUBJECT4" <FIX>
#	-"SUBJECT1 SUBJECT2 SUBJECT3 SUBJECT4" is list of subjects that can be processed
#	-<FIX>: if want to use fix to extract components then put in: classify  
#	-<FIX>: if want to classify and removal then put in: remove

# If nothing passed assume skip fix
var=${2:skip}

# Directory with subjects
cd /home/moreilly/Subjects

# Used for batch 
end=_ICA
mel=.ica
tmp=_struct_tmp
act=struct_recon-all
boldpathend='"'

# Loop through group of subjects
for SUBJECT in $1; do 

	# Put raw files in folder
	mkdir ./$SUBJECT/raw
	for rawfiles in `ls ./$SUBJECT`; do mv ./$SUBJECT/$rawfiles ./$SUBJECT/raw; done
	
	# Find run numbers 
	dcmunpack -src ./$SUBJECT/raw -scanonly ./$SUBJECT/scan_information.txt
	
	# If see nothing in bold folder- then might have to put these in manually 
	run_struct=`cat ./$SUBJECT/scan_information.txt | grep ' T1 ' | awk '{print $1}'`
	run_func=`cat ./$SUBJECT/scan_information.txt | grep RS_fMRI | awk '{print $1}'`

	# DICOM IMPORT
	dcmunpack -src ./$SUBJECT/raw/ -targ ./$SUBJECT \
		-run 2 $run_struct structural nii struct \
		-run 3 $run_func bold nii f
	
	# FREESURFER RECONSTRUCTION
	echo FREESURFER RECONSTRUCTION
	recon-all -i ./$SUBJECT/structural/*$run_struct/struct.nii -s $SUBJECT$tmp -all

	# Move reconstructed file
	mv ./$SUBJECT$tmp ./$SUBJECT/$act

	# SKULL STRIP
	bet ./$SUBJECT/structural/*$run_struct/struct.nii ./$SUBJECT/structural/*$run_struct/struct_brain.nii	

	# Edit design file with subject name 
	cp preprocess_design.fsf ./$SUBJECT
	
	# Find path to functionfsl
	boldpath=`find ./$SUBJECT/bold -name f.nii`
	sed -i -e "s/SUBJECT/$SUBJECT/" ./$SUBJECT/preprocess_design.fsf #replace names
	sed -i -e "s|/$SUBJECT/bold.*|${boldpath:1}$boldpathend|" ./$SUBJECT/preprocess_design.fsf #replace functional path
	sed -i -e "s/SUBJECT_ICA/$SUBJECT$end/" ./$SUBJECT/preprocess_design.fsf #replace ica path
	
	# MELODIC
	# motion correction, co-registration, high pass filtering, and ICA
	feat ./$SUBJECT/preprocess_design.fsf 
	
	# FIX 
        #if [ "$2" = classify ]; then 
	#	~/Documents/fix/fix -f ./$SUBJECT/$SUBJECT$end$mel
		# Change TRAIN.RData
	#	~/Documents/fix/fix -c ./$SUBJECT/$SUBJECT$end$mel TRAIN.RData 20
	#elif [ "$2" = remove ]; then
	#	~/Documents/fix/fix ./$SUBJECT/$SUBJECT$end$mel TRAIN.RData 20
	#fi

done 
