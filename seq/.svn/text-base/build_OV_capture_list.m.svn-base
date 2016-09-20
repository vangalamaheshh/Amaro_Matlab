function samples = build_OV_capture_list

fprintf('Getting list of capture samples.\n');
samples = {};
basedir = '/xchip/tcga_scratch/lawrence/ov';
d1 = dir([basedir '/*']);
for i=1:length(d1)
  if exist([basedir '/' d1(i).name '/capture'],'dir')
    shortname = d1(i).name;
    samples = [samples; 'ov/' shortname '/capture'];
  end
end
fprintf('Found %d samples\n',length(samples));
