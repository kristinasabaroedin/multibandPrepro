function []  = multiband_prepro(subject,runCrop,runSkullStrip,runFEATandMelodic,runFIX,runANTs,runSmoothing,runExtractTS)

    % Kristina Sabaroedin, Brain & Mental Health Laboratory, 2017
    
    clc;
    
    subprepro = ['/projects/kg98/kristina/GenofCog/datadir/derivatives/',subject];
    
    if exist(subprepro) == 0;
        sprintf('%s: Initialising derivatives folder\n', subject)
        mkdir(subprepro)
    end
    
    cd(subprepro)
    
    if exist('cfg.mat') == 2;
        load('cfg.mat')
        display('load cfg.mat')
    else
        display('Initialising cfg.mat')
    end
    
    cfg.Date = date;

    cfg.runCrop = runCrop;
    cfg.runSkullStrip = runSkullStrip;
    cfg.runFEATandMelodic = runFEATandMelodic;
    cfg.runFIX = runFIX;
    cfg.runANTs = runANTs;
    cfg.runSmoothing = runSmoothing;
    cfg.runExtractTS = runExtractTS;
    
    
    % -----------------------------------------------------------------------------------------------------------------------------
    % Before running FIX
    % Make sure that all required R libraries are installed (I think FIX
    % will prompt you to install the libraries, or it might fail and the
    % prompts will appear on the log file)
    % M3 users: please create an .Renviron in your home directory and replace
    % the path to the R package library that you installed
    % This ensures that R points to the libraries installed on your home
    % directory, instead of the library on 3 directory
    % Reference: https://csgillespie.github.io/efficientR/3-3-r-startup.html#renviron
    % All you need to do is create a text file, type the path to you
    % library (e.g.
    % R_LIBS=/home/ksabaroe/R/x86_64-pc-linux-gnu-library/3.3/partykit/libs)
    % and save it as .Renviron on your home directory
   
    % Note: FIX runs with Matlab/r2016a
    
    % -----------------------------------------------------------------------------------------------------------------------------
    % Paths to scripts required to run preprocessing
    % -----------------------------------------------------------------------------------------------------------------------------

    cfg.scriptdir = '/projects/kg98/kristina/GenofCog/scripts/prepro/';
    addpath(cfg.scriptdir)


    % -----------------------------------------------------------------------------------------------------------------------------
    % Paths to modules
    % -----------------------------------------------------------------------------------------------------------------------------

    % set FSL environment 
    cfg.fsldir = '/usr/local/fsl/5.0.9/fsl/bin/';
    setenv('FSLDIR',cfg.fsldir(1:end-4));
    setenv('FSLOUTPUTTYPE','NIFTI_GZ');
    setenv('LD_LIBRARY_PATH',[getenv('PATH'),getenv('LD_LIBRARY_PATH'),':/usr/lib/fsl/5.0'])

    % Directory of FIX
    cfg.fixdir = '/usr/local/fix/1.064/bin/fix/';
    setenv('FIXDIR',cfg.fixdir);
    
    
    % ANTs
    cfg.antsdir = '/usr/local/ants/2.2.0/bin/';
    setenv('ANTSPATH',cfg.antsdir);
	cfg.antsscriptsdir = '/usr/local/ants/2.2.0/bin/Scripts/';

    % Directory to AFNI functions
    cfg.afnidir = '/usr/local/afni/16.2.16/';
    addpath(cfg.afnidir)


    % -----------------------------------------------------------------------------------------------------------------------------
    % Paths to subject's data
    % -----------------------------------------------------------------------------------------------------------------------------
    
    cfg.subject = subject; 
    
    % where the raw data are
    cfg.rawdatadir = '/projects/kg98/kristina/GenofCog/datadir/rawdata/';

    % where the raw T1 is
    cfg.rawt1dir = [cfg.rawdatadir,subject,'/anat/'];

    % where the raw epi is
    cfg.rawepidir = [cfg.rawdatadir,subject,'/func/'];
    % Filename of raw epi
    cfg.rawepi = [subject,'_task-rest_bold.nii.gz'];

    % directory of derivatives data
    cfg.derivativesdir = '/projects/kg98/kristina/GenofCog/datadir/derivatives/';

    % Directory of preprocessed files
    % anatomical scan; this is where raw skull-stripped T1 is saved
    cfg.t1prepro = [cfg.derivativesdir, subject, '/anat/'];
    % general prepro; output from FEAT
    cfg.preprodir = [cfg.derivativesdir,subject,'/prepro/'];

    % filename of raw T1
    cfg.rawt1 = [subject,'_T1w.nii.gz'];
    % filename of cropped T1
    cfg.croppedt1 = [subject,'_T1w_crop.nii.gz']; 
    % filename of skull-stripped T1
    cfg.t1 = [subject,'_T1w_crop_brain.nii'];

    % Raw resting-state functional file is called through fsl’s design.fsf file, and will be processed using FEAT

    cfg.removeNoise = 'ICA-FIX';

    % -----------------------------------------------------------------------------------------------------------------------------
    % Set parameters
    % -----------------------------------------------------------------------------------------------------------------------------
    cfg.project = 'Genetics of Cognition';
    cfg.tN = 616;
    cfg.TR = 0.754;

    % -----------------------------------------------------------------------------------------------------------------------------
    % Set subject anatomical folder in derivatives directory 
    % -----------------------------------------------------------------------------------------------------------------------------
    if exist(cfg.t1prepro) == 0
        sprintf('%s: Initialising anatomical derivatives folder\n', subject)
        mkdir(cfg.t1prepro)
    end

    % -----------------------------------------------------------------------------------------------------------------------------
    % Set subject preprocessing folder in derivatives directory 
    % -----------------------------------------------------------------------------------------------------------------------------
    if exist(cfg.preprodir) == 0
        sprintf('%s: Initialising preprocessing folder\n', subject)
        mkdir(cfg.preprodir)
    end

    % -----------------------------------------------------------------------------------------------------------------------------
    % Crop out neck from raw T1
    % -----------------------------------------------------------------------------------------------------------------------------
    
    if runCrop == 1

        cd(cfg.rawt1dir)

        sprintf('%s: Cropping out neck from raw T1\n', subject)
        system([cfg.fsldir,'robustfov -i ',cfg.rawt1,' -r ', cfg.croppedt1]);
        movefile(cfg.croppedt1, cfg.t1prepro);
        display('Neck is cropped')

    end

    % -----------------------------------------------------------------------------------------------------------------------------
    % Skull-strip cropped T1
    % -----------------------------------------------------------------------------------------------------------------------------
    
    if runSkullStrip == 1

        cd(cfg.t1prepro)
        
        SkullStrip = {'ANTsBrainExtraction', 'BET'};
        WhichSkullStrip = SkullStrip{1};

		switch WhichSkullStrip
	
            case 'ANTsBrainExtraction'
                % Takes a longer time to run, requires a template and a mask,
                % but more robust than BET. Template and mask used here were
                % recommended by the creators of the software (for extracting healthy adult
                % brains)
                sprintf('%s: Performing ANTs Brain Extraction\n', subject)
                cfg.antsbrainextracttemplate = '/projects/kg98/kristina/templates/Oasis/T_template0.nii.gz';
                cfg.antsbrainextractmask = '/projects/kg98/kristina/templates/Oasis/T_template0_BrainCerebellumProbabilityMask.nii.gz';
                system([cfg.antsscriptsdir,'antsBrainExtraction.sh -d 3 -a ', cfg.croppedt1, ' -e ', cfg.antsbrainextracttemplate, ' -m', cfg.antsbrainextractmask, ' -o ANTs']);
                movefile('ANTsBrainExtractionBrain.nii.gz', [cfg.t1,'.gz']);
                movefile('ANTsBrainExtractionMask.nii.gz', [subject,'_T1w_crop_brain_mask.nii.gz']);
                display('T1 is skull stripped')

            case 'BET'
                % BET parameters used here was tailored to the Gen of Cog data,
                % check result to see if these parameters work for your data.
                % You might have to adjust the flags/values.
                sprintf('%s: Performing BET with -f 0.3 -m -R -B\n', subject)
                system([cfg.fsldir,'bet ',cfg.croppedt1,' crop_brain -f 0.3 -m -R -B']);
                movefile('crop_brain.nii.gz', [cfg.t1,'.gz']);
                movefile('crop_brain_mask.nii.gz', [subject,'_crop_brain_mask.nii.gz']); 
                display('T1 is BETted')
        end
    end

    % -----------------------------------------------------------------------------------------------------------------------------
    % Run FEAT and Melodic
    % Requires a set up design.fsf file that was manually saved from the FEAT gui
    % FEAT parameters
        % Delete first 4 volumes of epi
        % High pass filter: 75s
        % Smoothing: 3mm
        % FLIRT: 12 dof, BBR: input image is skull-stripped T1, reference is MNI 2mm T1 brain template
        % Non-linear registration turned ON, warping set to 10mm (default)
    % -----------------------------------------------------------------------------------------------------------------------------
   
    % Path to new FEAT directory (FEAT always creates a new folder every time it runs)
    cfg.featdir = [cfg.derivativesdir,subject,'/prepro.feat/'];
    
    % Main FEAT epi output
    cfg.featEpi = 'filtered_func_data.nii';
    
    meanEPIunsmoothed = 'mean_func_unsmoothed.nii.gz';
    
    if runFEATandMelodic == 1

        if exist(cfg.featdir) == 7;
            display('Reinitialising FEAT dir')
            rmdir(cfg.featdir, 's')
        else
            display('No existing FEAT directory - OK')
        end
            
        cd(cfg.derivativesdir)

        copyfile('design_master.fsf', cfg.preprodir)
        cd(cfg.preprodir)

        % This step below is superfluous, but I am pedantic about avoiding accidental alteration of the design_master.fsf file 
        movefile('design_master.fsf', 'design.fsf');

        % We want to edit design.fsf so the paths point to the current subject
        % create a variable that contains the content of the textfile
        buffer = fileread('design.fsf');
        % replaces ‘sub-000’ with subject within the buffer variable
        buffer = regexprep(buffer, 'sub-000', subject);
        % create and write a new textfile called design_run.fsf
        fid = fopen('design_run.fsf', 'w');
        % write the edited buffer content into new file
        fwrite(fid, buffer);
        % close text file
        fclose(fid);

        delete('design.fsf')
        
        % Run FEAT
        sprintf('%s: Running FEAT\n', subject)
        
        % For some reason FEAT only works on matlab when it's run as a
        % shell command
        system('sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/feat design_run.fsf"') 

        delete([cfg.preprodir,'design_run.fsf']);

    
        % -----------------------------------------------------------------------------------------------------------------------------
        % Retrieve unsmoothed mean EPI for ANTs normalisation
        % -----------------------------------------------------------------------------------------------------------------------------

        cd(cfg.featdir)

        sprintf('%s: Preparing unsmoothed mean epi for normalisation in ANTs\n', subject)

        copyfile([cfg.rawepidir,cfg.rawepi], cfg.featdir);

        system([cfg.fsldir, 'fslmaths ', cfg.rawepi, ' prefiltered_func_data -odt float'])

        display('Deleting 4 volumes')
        system([cfg.fsldir, 'fslroi prefiltered_func_data prefiltered_func_data 4 616'])

        display('Realigning data to middle volume')
        system([cfg.fsldir,'mcflirt -in prefiltered_func_data -out prefiltered_func_data_mcf -reffile example_func -rmsrel -rmsabs -spline_final'])


        display('Retrieving temporal mean of epi data')    
        system([cfg.fsldir, 'fslmaths prefiltered_func_data_mcf -Tmean mean_func_unsmoothed'])

        system([cfg.fsldir,'fslmaths mean_func_unsmoothed -mas mask mean_func_unsmoothed'])

        delete('prefiltered_func_data.nii.gz') 
        
        delete(cfg.rawepi) 

        rmdir('prefiltered_func_data_mcf.mat')
        
        delete('prefiltered_func_data_mcf*')

        display('Unsmoothed mean EPI is retrieved') 
    
        % -----------------------------------------------------------------------------------------------------------------------------
        % Run ICA Melodic 
        % -----------------------------------------------------------------------------------------------------------------------------

        % Directory of melodic output
        cfg.melodicdir = [cfg.featdir,'/filtered_func_data.ica'];

        cd(cfg.featdir)
        sprintf('%s: Running ICA-Melodic\n', subject)

        system('sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/melodic -i filtered_func_data.nii -o filtered_func_data.ica --nobet -m mask.nii --bgthreshold=5 --tr=0.754 --Ostats --report --mmthresh=0.5"') 
        
    end

    % -----------------------------------------------------------------------------------------------------------------------------
    % Run ICA FIX
    % Make sure training data already exists 
    % -----------------------------------------------------------------------------------------------------------------------------
    
    if runFIX == 1

        cfg.melodicdir = [cfg.featdir,'/filtered_func_data.ica'];
        
        cfg.trainingdata = '/projects/kg98/kristina/GenofCog/training/Training.RData';

        sprintf('%s: Running ICA-FIX\n', subject)
        
        cd(cfg.scriptdir);
       
        % FIX is called using an external shell script. This is because I
        % can't load R into matlab, which is required for FIX.
        
        % We want to edit the shell script so the paths point to the current subject
        % create a variable that contains the content of the textfile
        buffer = fileread('fix.sh');
        % replaces ‘sub-000’ with subject within the buffer variable
        buffer = regexprep(buffer, 'sub-000', subject);
        % create and write a new textfile called fix_run.sh
        fid = fopen(['fix_run_',subject,'.sh'], 'w');
        % write the edited buffer content into new file
        fwrite(fid, buffer);
        % close text file
        fclose(fid);
       
        system(['sh -c ". $bash /projects/kg98/kristina/GenofCog/scripts/prepro/fix_run_',subject,'.sh "'])

		delete(['fix_run_',subject,'.sh']);
    end
    
    % FIX epi output
    cfg.fixEpi = 'filtered_func_data_clean.nii.gz';
    % -----------------------------------------------------------------------------------------------------------------------------
    % Run ANTs registration
       
    % -----------------------------------------------------------------------------------------------------------------------------
    
    % Path to where the normalised images will be stored
    cfg.regdir = [cfg.featdir,'ants/'] 
    cfg.preprocesseddir = [cfg.featdir,'preprocessed/']
   
	cfg.meanEPIunsmoothed = meanEPIunsmoothed;

    if runANTs == 1

        display('Running ANTs')
        
        % Set parameters for ANTs registration
        % EPI input is cfg.fixEpi
        % T1 input is cfg.t1
        % Single EPI volume to be used in EPI to T1 registration is meanEPIunsmoothed

        % Path to the directory where segmented tissue files are stored
        % These tissue files were retrieved during FIX, in which FSL FAST segmented the T1 image into three tissue types
        fastdir = [cfg.featdir,'fix/'];

        % Tissue types
        cfg.gm = 'fastsg_pve_1.nii.gz';
        cfg.wm = 'fastsg_pve_2.nii.gz';
       	cfg.csf = 'fastsg_pve_0.nii.gz';

        % MNI template to be used in registration
        cfg.mni_template = '/projects/kg98/kristina/templates/MNI152_T1_2mm_brain.nii';

        cd(cfg.featdir)
        mkdir (cfg.regdir);
        copyfile(cfg.fixEpi, cfg.regdir)

        cd(cfg.regdir);
        
        sprintf('%s: Running ANTs registration\n', subject)

        SpatialNormalisationANTs([cfg.regdir, cfg.fixEpi],[cfg.featdir,cfg.meanEPIunsmoothed],...          
                [cfg.t1prepro,cfg.t1,'.gz'],...
                [fastdir,cfg.gm],...
                [fastdir,cfg.wm],...
                [fastdir,cfg.csf],...
                cfg.mni_template,cfg.antsdir,cfg.antsscriptsdir)

        % delete 
        delete(cfg.fixEpi);
        
        % rename files
        disp('renaming outputs');
        movefile (['w',cfg.fixEpi],[subject,'_',cfg.fixEpi(1:end-7),'_mni.nii.gz'])

        movefile (['w',cfg.gm], [subject,'_tissue_gm_mni.nii.gz'])
        movefile (['w',cfg.wm], [subject,'_tissue_wm_mni.nii.gz'])
        movefile (['w',cfg.csf], [subject,'_tissue_csf_mni.nii.gz'])

        movefile (['w',cfg.t1,'.gz'], [cfg.t1(1:end-4),'_mni.nii.gz'])
        
        % move files
        disp('Moving files to preprocessed directory')
        
        mkdir(cfg.preprocesseddir)
        
		% move epi
		movefile([cfg.regdir, subject,'_',cfg.fixEpi(1:end-7),'_mni.nii.gz'], cfg.preprocesseddir)
		% move t1
        movefile([cfg.regdir, cfg.t1(1:end-4),'_mni.nii.gz'], cfg.preprocesseddir)
		% move tissue files
        movefile([cfg.regdir,subject,'_tissue_gm_mni.nii.gz'], cfg.preprocesseddir)
        movefile([cfg.regdir,subject,'_tissue_wm_mni.nii.gz'], cfg.preprocesseddir)
        movefile([cfg.regdir,subject,'_tissue_csf_mni.nii.gz'], cfg.preprocesseddir)
        
        display('done')
     
    end

    cfg.normEpi = [subject,'_filtered_func_data_clean_mni.nii.gz']; 
	cfg.normT1 = [cfg.t1(1:end-4),'_mni.nii.gz'];
	cfg.normGM = [subject,'_tissue_gm_mni.nii.gz'];
	cfg.normWM = [subject,'_tissue_wm_mni.nii.gz'];
	cfg.normCSF = [subject,'_tissue_csf_mni.nii.gz'];
    
    % -----------------------------------------------------------------------------------------------------------------------------
    % Smooth normalised EPI
    % -----------------------------------------------------------------------------------------------------------------------------

    if runSmoothing == 1

        cd(cfg.preprocesseddir)
        
        % gunzip normEpi for AFNI smoothing
        %gunzip([cfg.normEpi,'.gz']);
        %delete([cfg.normEpi,'.gz']);

		cfg.smoothingkernel = 6;

        system([cfg.afnidir,'3dBlurToFWHM -FWHM ', num2str(cfg.smoothingkernel),' -mask /projects/kg98/kristina/templates/MNI152_T1_2mm_brain_mask_dil.nii -prefix smooth -input ', cfg.normEpi]);
        % Make sure -mask is pointing to the correct mask path, or you can
        % set mask to -automask

        display('Converting AFNI format to Nifti')

        system([cfg.afnidir,'3dAFNItoNIFTI smooth+*'])

        gzip('smooth.nii');
        
        delete('smooth.nii');
        
        display('Renaming output file')

		cfg.epiprepro = [cfg.normEpi(1:end-7),'_smooth.nii.gz'];

        movefile('smooth.nii.gz', cfg.epiprepro)
        
        delete('smooth+*'); 
    
        
    elseif runSmoothing == 0
        
        cfg.epiprepro = [cfg.normEpi(1:end-7),'_smooth.nii.gz'];
        
    end


    % -----------------------------------------------------------------------------------------------------------------------------
    % Extract time-series
    % -----------------------------------------------------------------------------------------------------------------------------

    if runExtractTS == 1

        cd(cfg.preprocesseddir)

        % Parcellation file for time series extraction
        cfg.parcFiles = {'/projects/kg98/kristina/ROIs/Gordon/Gordon_MNI_222.nii',...
                        '/projects/kg98/kristina/ROIs/Power/Power.nii',...
                        '/projects/kg98/kristina/ROIs/TriStri/TriStri.nii',...
                        '/projects/kg98/kristina/ROIs/ROIspheres/DiMartino.nii'};

        cfg.parcWeightGM = {'yes',...
                            'yes',...
                            'no',...
                            'no'};

        % Set input image for time series extraction
        cfg.ExtractIn = cfg.epiprepro;


        % Initialise roi time series variable
        cfg.roiTS = [];

            % Loop over parcellation files
            for i = 1:length(cfg.parcFiles)
                % Set parcellation file
                cfg.parcFile = cfg.parcFiles{i};
                cfg.parcName = cfg.parcFile((find(cfg.parcFile =='/',1,'last')+1:end-4));
                % set GM weight
                cfg.weightGM = cfg.parcWeightGM{i};
                % extract time series 
                cfg.roiTS{i} = prepro_extractTS_FSL(cfg);
            end

            cfg = rmfield(cfg,{'parcFile','weightGM','parcName','gm'});
        
                
            
           
        end

    % Save defined parameters and variables in a matlab file    
   
    cd([cfg.derivativesdir,subject])
   
    save('cfg.mat', 'cfg')

	rmdir('prepro')
   
    fprintf(1, '\t\t  %s: Preprocessing complete! \n', cfg.subject)

end


