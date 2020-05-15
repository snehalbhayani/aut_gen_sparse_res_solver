function [monsinstance] = subsmon(varvalues)
extendedbasis=[1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 6;1 1 1 2 2 3 0 0 0 0 0 0 0 1 1 1 1 1 2 2 2 3 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 3 0 0 0 0 0 0 1 1 1 2 0 0 0 1 0;0 1 2 0 1 0 0 0 1 1 2 2 3 0 0 1 1 2 0 0 1 0 0 0 0 1 1 1 2 2 3 0 0 0 1 1 2 0 0 1 0 0 0 0 1 1 2 0 0 1 0 0 0 1 0 0;2 1 0 1 0 0 2 3 1 2 0 1 0 1 2 0 1 0 0 1 0 0 1 2 3 0 1 2 0 1 0 0 1 2 0 1 0 0 1 0 0 0 1 2 0 1 0 0 1 0 0 0 1 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
for r=1:size(extendedbasis,2) 
temp=1; 
for j=1:size(extendedbasis,1) 
temp = temp * varvalues(j)^extendedbasis(j,r); 
end 
monsinstance(r,1) = temp; 
end
end
