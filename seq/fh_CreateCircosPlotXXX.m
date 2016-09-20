function fh_CreateCircosPlotXXX(libdir,circos_executable,perl_lib,individual_name,dRanger_file,karyotype_file,seg_file,refdir,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'need_to_chmod',true);
P = impose_default_value(P,'output_format','png');  % can be png or svg
P = impose_default_value(P,'confdir',libdir);
if exist('karyotype_file','var')
	% check karyotype_file for full path or just filename - if just a filename then use libdir
	[pathstr, name, ext] = fileparts(karyotype_file);
	if (length(pathstr)<1)
		P = impose_default_value(P,'karyotype_file',[libdir '/' karyotype_file]); 
	else
	 	P = impose_default_value(P,'karyotype_file',karyotype_file); 
	end
end

if exist('refdir','var')
	% reference info file
	P.refdir = refdir;
    ReferenceInfoObj.init(P.refdir)
    P.build=ReferenceInfoObj.getBuild())
end

        
P.dRanger_results_filespec = dRanger_file;
if exist('seg_file','var'), P.segfile_filespec = seg_file; end
P.CIRCOS_outname = [individual_name '.' P.output_format];
P.executable_name = circos_executable;
P.perl_modules_library = perl_lib;

if P.need_to_chmod
  cmd = ['chmod +x ' circos_executable];
  fprintf('Executing command: %s\n',cmd);
  try
    system(cmd);
  catch me
    fprintf('chmod failed\n');
  end
end

P = impose_default_value(P,'filter_by_somatic_score',true);
P = impose_default_value(P,'somatic_score_cutoff',4);

dRanger_draw_circos(individual_name,P);

if ~exist(P.CIRCOS_outname,'file')
  fprintf('CreateCircosPlot failed\n');
  cmd = ['cp ' libdir 'failure.png ./' P.CIRCOS_outname];  
  system(cmd);
end


