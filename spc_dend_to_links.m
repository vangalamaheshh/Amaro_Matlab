function lnk=spc_dend_to_links(d)
%SPC_DEND_TO_LINKS creates a links matrix from a dend matrix
%	LNK=SPC_DEND_TO_LINKS(D) creates a matrix of links describing
%	the dendrogram. Each link represents a binary split in the tree.
%	The first column in a level index that is used to represent a 
%	multi-son node. The links matrix has 5 columns : 
% 	   Left-start Left-end Right-start Right-end Level
%	and it start from heighest level i.e. bottom of tree
%
%	See also SPC_DEND
%

T=size(d,2)-2;
N=size(d,1);
lnk=[];
for i=1:T
%	[ i size(lnk,1) ]
	fs = 1;
	fe = 1;
	while ( fs <= N )
  		fcur = d(fs,i);
		while ( d(fe,i) == fcur ) & fe < N
			fe = fe+1;
		end
	
		if ( d(fe,i) ~= fcur )
			fe = fe-1;
		end
%		[ fs fe ]
%		pause(-1);
		% the father is from fs to fe
		% go over all father and find sons
		ss=fs;
		se=fs;
		while ( ss <= fe )
			scur = d(ss,i+1);
			while (d(se,i+1) == scur) & se < fe
				se = se +1;
			end
			if ( d(se,i+1) ~= scur )
				se = se -1;
			end
			% son is from ss to se
			if ( se ~= fe )
%				[ i ss se se+1 fe ]
				lnk = [lnk; ss se se+1 fe (i-1) ];
			end
			ss = se +1;
			se = ss;
		end		
		fs = fe+1;
		fe = fs;
	end;
end;	

lnk=flipdim(lnk,1);

