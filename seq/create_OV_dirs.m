function create_OV_dirs(X,subdirname)

require_fields(X,{'patient','bamname','bainame','isnorm'});
if ~exist('subdirname','var'), error('Must supply subdirname'); end

  function assertdir(dirname)
    existed = exist(dirname,'dir');
    if ~existed
      all_done = false;
      fprintf('Creating directory %s\n',dirname);
      success = mkdir(dirname);
      if ~success, error('Error creating directory %s!',dirname); end
    end
  end

  function assertlink(src,dest);
    existed = false;
    if exist(dest)
      d = dir(dest);
      s = dir(src);
      if s.bytes==d.bytes && s.datenum==d.datenum
        existed = true;
      else
        fprintf('Deleting out-of-date link %s\n',dest);
        result = system(['rm ' dest]);
        if result~=0, error('Error deleting %s!',dest); end
      end
    else
      % check to see if a "dead link" exists
      result = ls(dest);
      if contains(result,dest)
        fprintf('Deleting dead link %s\n',dest);
        result = system(['rm ' dest]);
        if result~=0, error('Error deleting %s!',dest); end
      end
    end
    if ~existed
      all_done = false;
      fprintf('Creating link %s\n',dest);
      result = system(['ln -s ' src ' ' dest]);
      if result~=0, error('Error creating link %s!',dest); end
    end
  end

mikedir = '/xchip/tcga_scratch/lawrence/ov';
ngdir = '/xchip/tcga_scratch/ng';
all_done = true;
p = unique(X.patient);
for i=1:length(p)
  tidx = find(strcmp(X.patient,p{i}) & ~X.isnorm);
  nidx = find(strcmp(X.patient,p{i}) & X.isnorm);
  if isempty(tidx) | isempty(nidx)
    fprintf('Skipping %s, which does not have tumor and normal bam\n',p{i});
  else
    % Mike's dirs
    destdir = [mikedir '/' p{i} '/' subdirname];
    assertdir(destdir);
    assertlink(X.bamname{tidx}, [destdir '/tumor.bam']);
    assertlink(X.bamname{nidx}, [destdir '/normal.bam']);
    assertlink(X.bainame{tidx}, [destdir '/tumor.bam.bai']);
    assertlink(X.bainame{nidx}, [destdir '/normal.bam.bai']);

    % ng dir structure
    destdir = [ngdir '/OV-' p{i}];
    assertdir(destdir);
    capdir = [destdir '/' subdirname];
    assertdir(capdir);
    bamdir = [capdir '/bam'];
    assertdir(bamdir);
    assertlink(X.bamname{tidx}, [bamdir '/tumor.bam']);
    assertlink(X.bamname{nidx}, [bamdir '/normal.bam']);
    assertlink(X.bainame{tidx}, [bamdir '/tumor.bam.bai']);
    assertlink(X.bainame{nidx}, [bamdir '/normal.bam.bai']);
  end
end

if all_done
  fprintf('All directories and links are already up-to-date.\n');
end


end % main function
