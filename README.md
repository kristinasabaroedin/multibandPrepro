# multibandPrepro
preprocessing scripts for multiband resting-state fmri
Main script: multiband_prepro.m

	This script is made to work with BIDS directory and file naming structure

	Raw unprocess nifti files should be under a rawdata folder (i.e. /projects/kg98/username/yourprojectname/rawdata/)

	All processed and intermediate files should go to a derivatives folder (i.e. /projects/kg98/username/yourprojectname/derivatives/)


	To test, this script can be run directly from matlab.

	Example (Matlab) if running all preprocessing steps: multiband_prepro('sub-987',1,1,1,1,1,1,1)

	Example slurm script: multiband_slurm.sh


FEAT

	FEAT uses a .fsf file to store its parameters

	.fsf file in this pipeline is named "design_master.fsf"

	Load this .fsf file into FSL FEAT to see what the parameters are, and change it to suit your data

	Leave the subject number to sub-000; the main scrpt will edit this number to match the subject ID that is being processed

	"design_master.fsf" should be stored in your derivatives directory. If you want to save it in a different directory, please make sure the multiband_prepro.m script is pointing to the correct path.

	A log of the FEAT processing will be in a "report_log.html" in your prepro.feat folder

	FEAT recreates a new .FEAT folder with "+" appended to it each time you re-run it. I.e. if you're rerunning FEAT on the same subject for the second time, the results will be saved in a folder called "prepro.feat+". Delete the existing .feat folder if you want to rerun FEAT.


ICA-FIX

	Training data: Training.RData

	Please edit directory to this training data in the multiband_prepro.m

	Shell script for running FIX: fix.sh

	Any error message from FIX will be stored in your script directory as a .txt file (usually named "errorLog.txt")


Normalisation to MNI

	Uses SpatialNormalisationANTs.m and antsRegistrationSyN.sh

	Transformation is done straight on the 4D data, which takes a lot of memory. If you are running this on Matlab using the MASSIVE desktop, operation might get killed. To make sure you don't run put of memory, make sure that you are connected to the LARGE desktop. This is not a problem if you are using slurm to run the script.


Extract TS

	Uses prepro_extractTS_FSL.m

	read.m is used in script prepro_extractTS_FSL.m to read parcellation files


Outputs

	Output and intermediate files will be in prepro.feat folder 

	Final preprocessed images should be in the preprocessed subfolder

	cfg.mat will be in the subject's parent derivative folder



Kristina Sabaroedin, 2017
