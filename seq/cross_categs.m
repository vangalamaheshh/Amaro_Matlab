function cross_categs(cset,varargin)
% cross_categs(cset,in_1,in_2,...,in_n,out)
%
% e.g. cross_categs(1:24,strand,zone,context65,allcateg)
%
% OR:
%
% cross_categs(build,cset,in_1,in_2,...,in_n,out)

beep off
if ~isdeployed
addpath /xchip/cga2/lawrence/cga/trunk/matlab/
addpath /xchip/cga2/lawrence/cga/trunk/matlab/seq
addpath /xchip/cga2/lawrence/cga/trunk/matlab/mike
addpath /xchip/cga2/lawrence/cga/trunk/matlab/seq/MutSig_2_dev

javaclasspath({...
'/xchip/cga2/lawrence/cga/trunk/matlab/seq/sam.jar';...
'/xchip/cga2/lawrence/cga/trunk/analysis_pipeline/tools/classes';...
});
import org.broadinstitute.cga.tools.seq.*;
import net.sf.samtools.*;
import java.io.*;
import java.lang.*;
end

build = '';
%keyboard
if ~isnumeric(cset)
  if strcmp(cset,'hg18')
    basedir = '/xchip/cga1/lawrence/db';
    build = 'hg18';
  elseif strcmp(cset,'hg19')
    basedir = '/xchip/cga1/lawrence/db/hg19';
    build = 'hg19'; 
  else 
    error('first parameter should be cset or build');
  end
cset = varargin{1};
%  varargin(1) = [];
%cset=1:24
else
  fprintf('Assuming hg18\n');
  build = 'hg18';
  basedir = '/xchip/cga1/lawrence/db';
end

categs = varargin;
if length(categs)<3, error('usage: cross_categs(cset,in_1,in_2,...,in_n,out)'); end
in_categs = categs(2:end-1);
out_categ = categs{end};

chrlen = load_chrlen(build);

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
c = cset(csetidx);
fprintf('\nchr%d\n',c);

%keyboard
if c<1 || c>24, error('c must be 1-24'); end

  len = chrlen(c);
  categ = zeros(len,1);
  categ_names = [];
  for i=1:length(in_categs)
    d = [basedir '/' in_categs{i}];
    in = load([d '/chr' num2str(c) '.mat']);
    f = fieldnames(in);
    if length(f)==1
      idx=1;
    else
      idx=find(strcmp(f,'categ'));
      if isempty(idx)
        fprintf('I don''t know which of these fields to use!');
        disp(f);
        fprintf('Please help me!\n');
        keyboard
      end
    end
    in = getfield(in,f{idx});
    in_categ_list = load_struct([d '/categs.txt'],'%f%s');
    in = 1+in-min(in_categ_list.num);

    if i==1
      categ = in;
      categ_names = in_categ_list.name;
    else
     if length(in)<length(categ), categ=categ(1:length(in)); end
     if length(categ)<length(in), in=in(1:length(categ)); end
     categ = ((categ-1)*slength(in_categ_list))+in;
      tmp = {};
      for j=1:length(categ_names)
        tmp = [tmp; regexprep(in_categ_list.name,'(.*)',[categ_names{j} ':$1'])];
      end
      categ_names = tmp;
    end
  end

  d = [basedir '/' out_categ];
  if ~exist(d,'dir'), mkdir(d); end

  if c==1
    tmp = [];
    tmp.num = (1:length(categ_names))';
    tmp.name = categ_names;
    save_struct(tmp,[d '/categs.txt']);
  end

  fname = [d '/chr' num2str(c) '.mat'];
  save(fname,'categ');
  f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
  fprintf(f,'%d\n',categ);
  fclose(f);

end
