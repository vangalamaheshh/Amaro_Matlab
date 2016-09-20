function check_cnqc(samples,params)
% check_cnqc(samples,params)
%
% non-firehose version of fh_CopyNumberQCReport
%
% Mike Lawrence 2009-10-29

if ~iscell(samples), samples={samples}; end
if ~exist('params','var'), params = []; end
params = impose_default_value(params,'auto',false);
params = impose_default_value(params,'flag1',2);

for i=1:length(samples), sample = samples{i};
  fprintf('Sample %d/%d: %s\n',i,length(samples),sample);
  individual_name = upper(regexprep(sample,'/','-'));
  direc = ['/xchip/tcga_scratch/lawrence/' sample];

%  try
    % load and process data
    [Lt Ln R T N Z] = load_cnqc(sample);

    % visualize
    draw_CNQC_plot(Lt,Ln,R,T,N,Z,individual_name);
    print('-dpng','-r100',[direc '/' individual_name '_CopyNumberQC.png']);

    % QC report: textfile
    L = generate_CNQC_report(Lt,Ln,individual_name);
    save_struct(L,[direc '/' individual_name '_CopyNumberQC_report.txt']);

    % QC report: HTML
    H = convert_CNQC_report_to_html(L,individual_name);
    save_textfile(H,[direc '/' individual_name '_CopyNumberQC_report.html']);

%  catch me, fprintf('FAILED: %s\n', individual_name); disp(me.message); disp(me); end
  if (params.auto), continue; end

  % show report
  fprintf('\n');for i=1:slength(L);fprintf('lane %d\t%s\t(%d)\t%s\t%s\t(%0.2f)\n',...
    i,L.tn{i},L.nreads(i),L.baitset{i},L.judgement{i},L.tumorseg_corr(i));end
  n_mixups = length(grep('MIXUP|CONTAM',L.judgement,1));
  fprintf('\nNumber of mixups: %.0f\n\n',n_mixups);

  if i<length(samples)
    fprintf('Type "return" to continue.\n');
    keyboard
  end
end

fprintf('Type "return" to end.\n');
keyboard

