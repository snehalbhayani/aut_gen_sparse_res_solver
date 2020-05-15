function [monsinstance] = subsmon(varvalues)
extendedbasis=[0 0 0 1 1 1 1 2 2;0 1 2 0 0 1 1 0 0;0 0 1 0 1 0 1 0 1;0 0 0 0 0 0 0 0 0]; 
for r=1:size(extendedbasis,2) 
temp=1; 
for j=1:size(extendedbasis,1) 
temp = temp * varvalues(j)^extendedbasis(j,r); 
end 
monsinstance(r,1) = temp; 
end
end
