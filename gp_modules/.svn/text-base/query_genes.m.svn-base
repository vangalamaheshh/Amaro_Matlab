function query_genes(a)
% query_genes  -list genelist with header -ui path to input struct -stats stats -regs regs -o output txt file name 

%a=handle_args({'list','ui','stats','regs', 'o'},varargin);
CL21=a.ui;
x=load(CL21)
fn=fieldnames(x);
D=getfield(x,fn{1});
regs=a.regs;
stats=a.stats;
load(regs);
load(stats);
load('~gadgetz/projects/snp/data/Refgene/hg17/ucsc_20070131/hg17_20070131.mat');
rg=add_chrn(rg);

genelist=a.list;
genes=read_dlm_file(genelist);

notfound={};
nomatch={};
wide_genes={};
narrow_genes={};
type={};

for i=2:length(genes)
    w=grep(strcat('^', genes{i}(1), '$'),{rg.symb},1);
    if ~isempty(w)
        % w is index, rg(w(1)) is the refseq entry
        % widx(i)=w(1);

        snps=find_snps(D,rg(w(1)).chrn,rg(w(1)).start,rg(w(1)).end,1);
        for k=1:2
            for j=1:length(regs{k})

                if regs{k}(j).chrn==rg(w(1)).chrn && regs{k}(j).peak_st < min(snps) && regs{k}(j).peak_en > max(snps)
                    narrow_genes(i)=genes{i}(1);
                    if k==1
                        type{i}='amplification';
                    else
                        type{i}='deletion';
                    end
                elseif regs{k}(j).chrn==rg(w(1)).chrn && regs{k}(j).peak_wide_st < min(snps) && regs{k}(j).peak_wide_en > max(snps)
                    wide_genes(i)=genes{i}(1);
                   if k==1
                        type{i}='amplification';
                    else
                        type{i}='deletion';
                    end
                else
                    notfound(i)=genes{i}(1);
                end
            end
        end
    else
        nomatch(i)=genes{i}(1);

    end
end

%keyboard
save('outfile.mat', 'genes', 'wide_genes', 'notfound', 'nomatch', 'narrow_genes') 

if isempty(a.o)
   error('Missing output file name');
end
outfile=a.o;
f=fopen(outfile,'w');

fprintf(f,['genes in wide region:' '\n']);
for i=1:length(wide_genes)
    if ~isempty(wide_genes{i})
fprintf(f,[wide_genes{i} '--'    type{i} '\n']);
    end
end

fprintf(f, ['genes in narrow region:' '\n']);
for i=1:length(narrow_genes)
    if ~isempty(narrow_genes{i})
    fprintf(f, [narrow_genes{i} '--'    type{i} '\n']);
    end
end

%disp('not found:');
%for i=1:length(notfound)
%disp(notfound{i});
%end

fprintf(f,  ['no match:' '\n']);
for i=1:length(nomatch)
fprintf(f, [nomatch{i} '\n']);
end
fclose(f);
