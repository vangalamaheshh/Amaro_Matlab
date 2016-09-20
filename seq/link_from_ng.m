function link_from_ng(samplemask,filemask,ngsubdir)
% function link_from_ng(samplemask,filemask,ngsubdir)
%
% e.g.
% link_from_ng('*/*/wgs','*zone*txt','cov')
%
% for each directory /xchip/tcga_scratch/lawrence/*/*/wgs
%   for each file *zone*txt in that directory
%     creates a link to that file from /xchip/tcga_scratch/ng/*-*/wgs/cov
%
% notes
%  (1)  */*/wgs will be automatically converted to *-*/wgs
%       e.g. ov/0751/wgs -> OV-0751/wgs
%
%  (2)  off-target links will be deleted and refreshed
%
% Mike Lawrence 2009-08-05

bd1 = '/xchip/tcga_scratch/lawrence';
bd2 = '/xchip/tcga_scratch/ng';

[status list] = system(['find ' bd1 '/' samplemask '/' filemask]);
if status ~=0, error('find failed'); end

list = text_to_lines(list);
for i=1:length(list)
  t = regexp(list{i},['^' bd1 '/(.*)/(.*)/(.*)/(.*)$'],'tokens');
  if isempty(t), error('token extract failed: %s\n',list{i}); end
  tumortype = t{1}{1};
  sample = t{1}{2};
  tech = t{1}{3};
  filename = t{1}{4};
  todir = [bd2 '/' upper(tumortype) '-' upper(sample) '/' tech '/' ngsubdir];
  if ~exist(todir), fprintf('Creating directory %s\n', todir); mkdir(todir); end
  link = [todir '/' filename];
  already_exists = false;
  if exist(link,'file')
    dto = dir(link);
    dfrom = dir(list{i});
    if ~strcmp(dto.name,dfrom.name) || dto.datenum~=dfrom.datenum
      fprintf('Deleting previous file %s\n', link);
      delete(link);
    else
      already_exists = true;
    end
  end
  if ~already_exists
    fprintf('Linking %s -> %s\n', list{i}, link);
    system(['ln -s ' list{i} ' ' link]);
  end
end



