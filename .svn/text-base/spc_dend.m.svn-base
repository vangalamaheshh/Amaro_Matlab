function dendlabs=spc_dend(labs)
%SPC_DEND orders labs in a dendogram fashion.
%	D=SPC_DEND(L) orders the labels and the points such that
%	points belonging to the same cluster are concecutive.  
%	This is done for all temperatures.
%
%	If there is a conflict, meaning two points which are in the same
%	cluster at a high temperature are NOT in the same cluster at a
%	lower temperature, the lower temperature wins, i.e a split can never
%	be connected at higher temp.
%
%	D has an additional three columns, the first column is all 1,
%	and the last two columns contain labels 1:N and the original index of the label.
%
%	See also SPC_READ_LABS,SPC_DRAW_DEND
% 


%load /users3/gaddyg/cluster/fmri/scr2/avg.con/all.labels

N=size(labs,1);
T=size(labs,2);
new = zeros(N,T+3);

% all data belong to first cluster at i=1
new(:,1)=ones(N,1);

new(:,2:T+1)=labs;
new(:,T+2)=(1:N)';
new(:,T+3)=(1:N)';
new2 = new;

% go over all t
for i=2:T+2,
	first = 1;
	last = 1;
%	disp('new  temprature');

% check for bad data
	% create matrix of labels for current transition
	nc1 = max(new(:,i-1)); % number of clusters in i-1
	nc2 = max(new(:,i)); % number of clusters in i
 
	q = sparse([],[],[],nc1,nc2);
	for j=1:N,
		q(new(j,i-1),new(j,i)) = q(new(j,i-1),new(j,i)) +1;
	end;


	% check that all new clusters had uniqe "fathers" 
%	tot = 0;
%	for j=1:nc2,
%		[dm,di]=max(q(:,j));
%		if ( sum(q(:,j)) ~= q(di,j) )
%			tot = tot + sum(q(:,j)) - max(q(di,j));
%			tmp = q(:,j);
%			tmp(di)=0;
%			j
%			dd = find( tmp> 0 )
%
%			% good for only 1 bad point
%			for k=1:N,
%				if ( new(k,i-1) == dd & new(k,i)== j )
%					new(k,size(new,2))
%				end;
%			end;
%		end;	
%	end;
%	if ( tot > 0 )
%		
%	end;
%	[ i nc1 tot ]
	
% go over all clusters in previous t
	j = 0;
	while ( first <= N )

		% find region of cluster j [first last]
		j = j + 1;
		cur = new(first,i-1);
		last = last +1;
		stop = 0;
		while (last <= N) & ( stop == 0 ),
			if ( new(last,i-1) == cur )
				last = last +1;
			else
				stop = 1;
			end;			
		end;
		last = last -1; 
		if ( last > first )
			[n,idx] = sort(new(first:last,i));
			new2(first:last,:)=new(idx+first-1,:);
		end;
		first = last+1;
		last = first;
	end;
%	[ j max(new(:,i-1)) ]
	new = new2;
end;

dendlabs=new2;


