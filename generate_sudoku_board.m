function [SB] = generate_sudoku_board( )
%generates some random numbers places them on the board and checks to see
%if its solvable. When it finds a solvable one. 

locs(:,1)=randi([1 9],1,12);
locs(:,2)=randi([1 9],1,12);

SB_matrix_test=zeros(9,9);
for i=1:size(locs,1)
    SB_matrix_test(locs(i,1),locs(i,2))=randi([1 9],1,1);
end

pass=0;

while ~pass
 for i=1:9
 rows=find(SB_matrix_test(i,:)~=0);
 
 end
end


end

