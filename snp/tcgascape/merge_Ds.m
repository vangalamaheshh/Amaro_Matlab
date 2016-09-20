function D = merge_Ds(D1,D2)
%merge_Ds - merge two SegArray-based D structs
%
% (very special purpose at this juncture)

% if either D struct is empty, nothing to merge, return the other
if isempty(D1)
    D = D2;
    return
end
if isempty(D2)
    D = D1;
    return
end

[~,idx1,idx2] = match_num_sets(1e11*D1.chrn+D1.pos,1e11*D2.chrn+D2.pos);
%! match_num_sets does not return array if args are unique sets! TODO: fix!!!
%! M = match_num_sets(1e11*D1.chrn+D1.pos,1e11*D2.chrn+D2.pos);
%! [idx1,idx2] = find(M);

% assert that sample identifiers are unique
[~,sx1,sx2] = match_string_sets_hash(D1.sdesc,D2.sdesc);
if ~isempty(sx1)
    error('%d identifiers are not unique!',length(sx2));
end

lx1 = false(size(D1.pos,1),1);
lx1(idx1) = true;
D1 = reorder_D_rows(D1,SegArray(lx1));

lx2 = false(size(D2.pos,1),1);
lx2(idx2) = true;
D2 = reorder_D_rows(D2,SegArray(lx2));

% create structure and merge fields
D = struct;
D.sdesc =    [D1.sdesc;D2.sdesc];
D.chrn =      D1.chrn;
D.pos =       D1.pos;
D.dat =      [D1.dat,D2.dat];
if isfield(D1,'sis') && isfield(D2,'sis')
    D.sis =      [D1.sis;D2.sis];
end
D.suppress_history = true;
D.islog = true;
D.isMB = false;
