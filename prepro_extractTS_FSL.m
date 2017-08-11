function roiTS = prepro_extractTS_FSL(cfg)
    % This function will extract resting state time series values from a parcellation image in MNI space.
    % It will also weight the averaging by grey matter probability
    % 
    % ------
    % INPUTS
    % ------
    % cfg.subject       - string containing subject ID = e.g., '1008.2.48.009'
    %
    % cfg.removeNoise   - string defining noise strategy used for preprocessing (see prepro_noise.m)
    % 
    % cfg.parcFile      - location and name of MNI parcellation nifti file. e.g., /path/to/dir/parcellation.nii
    % 
    % cfg.t1dir         - location of grey matter probability mask
    % cfg.normGM            - name of grey matter prob mask
    %
    % Linden Parkes, Brain & Mental Health Laboratory, 2016
    % ------------------------------------------------------------------------------

    fprintf('\n\t\t ----- Times series extraction ----- \n\n');
    fprintf(1, '\t\t Subject: %s \n', cfg.subject)
    fprintf(1, '\t\t Noise correction: %s \n', cfg.removeNoise);
    fprintf(1, '\t\t Input file: %s \n', cfg.ExtractIn);
    fprintf(1, '\t\t Parcellation file: %s \n',cfg.parcFile);
    fprintf(1, '\t\t Parcellation name: %s \n',cfg.parcName);
    fprintf(1, '\t\t GM weighting: %s \n', cfg.weightGM);
    fprintf('\n\t\t ----------------------------------- \n\n');

    % ------------------------------------------------------------------------------
    % ROIs
    % ------------------------------------------------------------------------------
    [parc_hdr,parc] = read(cfg.parcFile);
    numROIs = numel(unique(parc))-1;

    % ------------------------------------------------------------------------------
    % Script
    % ------------------------------------------------------------------------------
    cd(cfg.preprocesseddir)

    switch cfg.weightGM
        case 'yes'
            % multiply epi by gm prob
            system([cfg.fsldir,'fslmaths ',cfg.ExtractIn,' -mul ',cfg.preprocesseddir,cfg.normGM,' epi_w']);

            % take the average of the weighted epi data for each parcel
            system([cfg.fsldir,'fslmeants -i epi_w.nii* --label=',cfg.parcFile,' -o epi_parc_mean.txt']);
        
            % take the average of the gm probabilities for each parcel
            system([cfg.fsldir,'fslmeants -i ',cfg.preprocesseddir,cfg.normGM,' --label=',cfg.parcFile,' -o gm_parc_mean.txt']);

            x = dlmread('epi_parc_mean.txt');
            y = dlmread('gm_parc_mean.txt');

            for i = 1:numROIs
                roiTS(:,i) = x(:,i)/y(i);
            end

            delete('epi_w.nii*','epi_parc_mean.txt','gm_parc_mean.txt')
            
        case 'no'
            % take the average of the weighted epi data for each parcel
            system([cfg.fsldir,'fslmeants -i ',cfg.ExtractIn,' --label=',cfg.parcFile,' -o roiTS_',cfg.parcName,'.txt']);            
            roiTS = dlmread(['roiTS_',cfg.parcName,'.txt']);
            delete('roiTS_*.txt')
    end 
    

    fprintf('\n\t\t ----- ROI time series extracted ----- \n\n');
end
