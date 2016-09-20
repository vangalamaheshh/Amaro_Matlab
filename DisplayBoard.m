function [ ] = DisplayBoard( Board_Structure )
%Takes in the Board structure with the 9 quadrants and displays them as one
%unit for the user to manipulate;

Board_matrix(1:3,1:3)=Board_Structure.top1;
Board_matrix(1:3,4:6)=Board_Structure.top2;
Board_matrix(1:3,7:9)=Board_Structure.top3;

Board_matrix(4:6,1:3)=Board_Structure.mid1;
Board_matrix(4:6,4:6)=Board_Structure.mid2;
Board_matrix(4:6,7:9)=Board_Structure.mid3;

Board_matrix(7:9,1:3)=Board_Structure.bot1;
Board_matrix(7:9,4:6)=Board_Structure.bot2;
Board_matrix(7:9,7:9)=Board_Structure.bot3;

Board_Cell=cell(11,11);
for i=1:size(Board_matrix,1)+2
    for j=1:size(Board_matrix,2)+2
    if i==4 || j==4 || i==8 || j==8 
        Board_Cell{i,j}='--';
    if j==4 || j==8
        Board_Cell{i,j}='|';
    else    
    end
    end
    end
end
for i=1:3
    for j=1:3
   Board_Cell{i,j}=Board_matrix(i,j);
   Board_Cell{i,j+4}=Board_matrix(i,j+3);
   Board_Cell{i,j+8}=Board_matrix(i,j+6);
   
   Board_Cell{i,j}=Board_matrix(i,j);
   Board_Cell{i+4,j}=Board_matrix(i+3,j);
   Board_Cell{i+8,j}=Board_matrix(i+6,j);
   
   
   Board_Cell{i,j}=Board_matrix(i,j);
   Board_Cell{i+4,j+4}=Board_matrix(i+3,j+3);
   Board_Cell{i+8,j+8}=Board_matrix(i+6,j+6);
   
   Board_Cell{i+4,j+8}=Board_matrix(i+3,j+6);
   Board_Cell{i+8,j+4}=Board_matrix(i+6,j+3);
    end
end

%disp(Board_matrix)
disp(Board_Cell)
end

