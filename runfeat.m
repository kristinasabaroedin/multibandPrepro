fsl_path = '/usr/local/fsl/5.0.9/fsl/';
setenv('FSLDIR',fsl_path)
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
curpath = getenv('PATH');
setenv('PATH',sprintf('%s:%s',fullfile(fsl_path,'bin'),curpath));

system('sh -c ". ${FSLDIR}feat"')

For example:
system('sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/feat design_run.fsf"') 

system('sh -c ". /usr/local/fsl/5.0.9/etc/fslconf/fsl.sh;/usr/local/fsl/5.0.9/bin/feat design_run.fsf"')