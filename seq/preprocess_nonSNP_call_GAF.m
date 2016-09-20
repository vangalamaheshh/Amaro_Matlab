function I = preprocess_nonSNP_call(GAF, chr, st, en, change)
% preprocess_nonSNP_call(R, build, chr, st, en, change,P)
%
% Given Refseq database, find's overlapping transcripts, adjusts indel coordinates, and returns 
% necessary variables for find_indel_in_refseq.m
%
% Updated to handle DNPs
%
% Alex Ramos May 2010

	build = 'hg19';
	%fprintf('%d\t%d',st,en)
	
%	if ~exist('P','var'), P=[]; end
%	P = impose_default_value(P,'adjust_coordinates_for_Broad',false);
	
	
	if ~isfield(GAF,'gene_start'), GAF.gene_start = GAF.tx_start; end
	if ~isfield(GAF,'gene_end'), GAF.gene_end = GAF.tx_end; end
	pos1 = st;  %variable for indel start
	pos2 = en;  %variable for indel end
	I = [];

	
	if strncmp(change, '+', 1)  
		I.classification = 'INS';
		I.change = change(2:end);
		
%		if P.adjust_coordinates_for_Broad
%			if st == en
%				en = en + 1;  %adjust insertion coordinates
%			else
%				error('Insertion start and end coordinates expected to be equal');
%			end
%		end
		
		
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
		
%		if P.adjust_coordinates_for_Broad
%			st = st+1;
%		end
		
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
	
	pos_stridx = find(strcmp(GAF.strand,'+'));
	neg_stridx = find(strcmp(GAF.strand,'-'));
	%if st == 66455993
	%	error('LOL')
	%end
	%fprintf('%d\t%d',I.pos1,I.pos2)
  	if ~exist('idx','var') || isempty(idx)
  		if xor(isnumeric(GAF.chr),isnumeric(chr))
    		error('Query and GAF chromosome formats incompatible');
  		else
  			if isnumeric(GAF.chr) %!A - indel
	  			chridx = find(GAF.chr==chr);
  			else
	    		if xor(strncmp(chr,'chr',3), strncmp(GAF.chr{1},'chr',3))
	      			error('Query and RefSeq chromosome formats incompatible');
	    		else %!A - indel
	    			chridx = find(strcmp(GAF.chr,chr));
				end
			end
			
			pos_stridx = intersect(chridx, pos_stridx); 
			neg_stridx = intersect(chridx, neg_stridx);
  			% idx1 - for indels within a gene
  			% idx2 - for indels overlapping an IGR|gene_start boundary
  			% idx3 - for indels overlapping an gene_end|IGR boundary
			% idx4 - for indels encompassing a whole gene
    		pos_idx1 = pos_stridx(GAF.gene_start(pos_stridx) <= pos1 & GAF.gene_end(pos_stridx) >= pos2);
    		pos_idx2 = pos_stridx(GAF.gene_start(pos_stridx) >= pos1 & GAF.gene_start(pos_stridx) <= pos2);
    		pos_idx3 = pos_stridx(GAF.gene_end(pos_stridx) >= pos1 & GAF.gene_end(pos_stridx) <= pos2);
			pos_idx4 = pos_stridx(GAF.gene_start(pos_stridx) >= pos1 & GAF.gene_end(pos_stridx) <= pos2);
			pos_idx = [pos_idx1; pos_idx2; pos_idx3; pos_idx4];
			
    		neg_idx1 = neg_stridx(GAF.gene_start(neg_stridx) >= pos1 & GAF.gene_end(neg_stridx) <= pos2);
    		neg_idx2 = neg_stridx(GAF.gene_start(neg_stridx) >= pos1 & GAF.gene_start(neg_stridx) <= pos2);
    		neg_idx3 = neg_stridx(GAF.gene_end(neg_stridx) >= pos1 & GAF.gene_end(neg_stridx) <= pos2);
			neg_idx4 = neg_stridx(GAF.gene_start(neg_stridx) <= pos2 & GAF.gene_end(neg_stridx) >= pos1);
			neg_idx = [neg_idx1; neg_idx2; neg_idx3; neg_idx4];			
			
			idx = [pos_idx; neg_idx];
			idx = unique(idx);
			
    	end
	end
	I.idx = idx;
  
end  %end prepocess_indel