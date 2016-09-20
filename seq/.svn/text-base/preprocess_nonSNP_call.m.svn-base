function I = preprocess_nonSNP_call(R, build, chr, st, en, change,P)
% preprocess_nonSNP_call(R, build, chr, st, en, change,P)
%
% Given Refseq database, find's overlapping transcripts, adjusts indel coordinates, and returns 
% necessary variables for find_indel_in_refseq.m
%
% Updated to handle DNPs
%
% Alex Ramos May 2010


	%fprintf('%d\t%d',st,en)
	
	if ~exist('P','var'), P=[]; end
	P = impose_default_value(P,'adjust_coordinates_for_Broad',false);
	
	
	if ~isfield(R,'gene_start'), R.gene_start = R.tx_start; end
	if ~isfield(R,'gene_end'), R.gene_end = R.tx_end; end
	pos1 = st;  %variable for indel start
	pos2 = en;  %variable for indel end
	I = [];

	
	if strncmp(change, '+', 1)  
		I.classification = 'INS';
		I.change = change(2:end);
		
		if P.adjust_coordinates_for_Broad
			if st == en
				en = en + 1;  %adjust insertion coordinates
			else
				error('Insertion start and end coordinates expected to be equal');
			end
		end
		
		
		if iscell(chr)
			error('chr should not be a cell');
		elseif isnumeric(chr)
			l = chrlabel(chr);
			I.genomechange = ['g.' l{1} ':'...
			num2str(st) '_' num2str(en) 'ins' I.change];
		elseif strncmp(chr,'chr',3)
			I.genomechange = ['g.' chr ':'...
			num2str(st) '_' num2str(en) 'ins' I.change];
		else
			I.genomechange = ['g.chr' chr ':'...
			num2str(st) '_' num2str(en) 'ins' I.change];
		end
		
		
	elseif strncmp(change, '-', 1)
		I.classification = 'DEL';
		I.change = change(2:end);
		
		if P.adjust_coordinates_for_Broad
			st = st+1;
		end
		
		
		if length(I.change) == 1
			outl = [num2str(st) 'del' I.change];
		else
			outl = [num2str(st) '_' num2str(en) 'del' I.change]; 
		end
		
		if iscell(chr)
			error('chr should not be a cell');
		elseif isnumeric(chr)
			l = chrlabel(chr);
			I.genomechange = ['g.' l{1} ':' outl];
		elseif strncmp(chr,'chr',3)
			I.genomechange = ['g.' chr ':'  outl];
		else
			I.genomechange = ['g.chr' chr ':' outl];
		end
		
	elseif strncmp(change, '~',1)
		old_bases = upper(genome_region(chr,st,en,build));
		I.change = change(2:end);
			if length(I.change) == 2
				I.classification = 'DNP';
			elseif length(I.change) == 3
				I.classification = 'TNP';
			elseif length(I.change) > 3
				I.classification = 'ONP';
			else
				error('Incorrect change length for DNP!');
			end
			
		if iscell(chr)
			error('chr should not be a cell');
		elseif isnumeric(chr)
			l = chrlabel(chr);
			I.genomechange = ['g.' l{1} ':'...
			num2str(st) '_' num2str(en) old_bases '>' I.change];
		elseif strncmp(chr,'chr',3)
			I.genomechange = ['g.' chr ':'...
			num2str(st) '_' num2str(en) old_bases '>' I.change];
		else
			I.genomechange = ['g.chr' chr ':'...
			num2str(st) '_' num2str(en) old_bases '>' I.change];
		end
	
	else
		error('change variable for indel not in correct format');
	end	
	
	I.pos1 = st;
	I.pos2 = en;
	%fprintf('%d\t%d',I.pos1,I.pos2)
  	if ~exist('idx','var') || isempty(idx)
  		if xor(isnumeric(R.chr),isnumeric(chr))
    		error('Query and RefSeq chromosome formats incompatible');
  		elseif isnumeric(R.chr) %!A - indel
  			% idx1 - for indels within a gene
  			% idx2 - for indels overlapping an IGR|gene_start boundary
  			% idx3 - for indels overlapping an gene_end|IGR boundary
			% idx4 - for indels encompassing a whole gene
    		idx1 = find(R.chr==chr & R.gene_start <= pos1 & R.gene_end >= pos2);
    		idx2 = find(R.chr==chr & R.gene_start >= pos1 & R.gene_start <= pos2);
    		idx3 = find(R.chr==chr & R.gene_end >= pos1 & R.gene_end <= pos2);
			idx4 = find(R.chr==chr & R.gene_start >= pos1 & R.gene_end <= pos2);
    		idx = union(union(idx1,idx2), union(idx3,idx4));
  		else
    		if xor(strncmp(chr,'chr',3), strncmp(R.chr{1},'chr',3))
      			error('Query and RefSeq chromosome formats incompatible');
    		else %!A - indel
    			idx1 = find(strcmp(R.chr,chr) & R.gene_start <= pos1 & R.gene_end >= pos2);
    			idx2 = find(strcmp(R.chr,chr) & R.gene_start >= pos1 & R.gene_start <= pos2);
    			idx3 = find(strcmp(R.chr,chr) & R.gene_end >= pos1 & R.gene_end <= pos2);
				idx4 = find(strcmp(R.chr,chr) & R.gene_start >= pos1 & R.gene_end <= pos2);
    			idx = union(union(idx1,idx2), union(idx3,idx4));
    		end
  		end
	end
	I.idx = idx;
  
end  %end prepocess_indel