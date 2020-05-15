
afile = fopen(strcat(foldername, '/gen_elim_matrices.m'),'w');
fprintf(afile, '%s \n', 'function [A,B] = gen_elim_matrices(A,B)');


thresh = 1e-10;
nonzeroacols = find(sum(abs(A)) < thresh);
candcols = nonzeroacols;
ext = 4;
%% Version 2
while ~isempty(nonzeroacols)
    
    retries = 200;
    colsdone=[];
    subtractors=[];
    
    Btemp = B(:,find(sum(abs(A))<thresh));
    fprintf(afile, '%s \n', strcat('[q,r] = qr(B(:,', mat2str(find(sum(abs(A))<thresh)),'));'));
    fprintf(afile, '%s \n', 'A=q\A; B=q\B;');
    [q, r] = qr(Btemp);
    A = q\ A;
    B = q\ B;
    
    candcols = find(sum(abs(A)) < thresh); 
    fprintf(afile, '%s\n', strcat('A(:,', mat2str(candcols),') =[];'));
    fprintf(afile, '%s\n', strcat('B(:,', mat2str(candcols),') =[];'));
    fprintf(afile, '%s\n', strcat('A(', mat2str([1:size(candcols,2)]),',:) =[];'));
    fprintf(afile, '%s\n', strcat('B(', mat2str([1:size(candcols,2)]),',:) =[];'));
    extendedbasis(:,candcols) = [];
    A(:,candcols) = [];
    B(:,candcols) = [];
    A([1:size(candcols,2)], :) = [];
    B([1:size(candcols,2)], :) = [];
    nonzeroacols = find(sum(abs(A)) < thresh);
end
% %% Version 1
% while ~isempty(nonzeroacols)
%     
%     retries = 200;
%     colsdone=[];
%     subtractors=[];
%     
%     Btemp = B(:,find(sum(abs(A))<thresh));
%     while ~isempty(candcols)
% %         break;
%         Btemp = B(:,find(sum(abs(A))<thresh));
%         imagesc(flipud(abs(Btemp)>thresh));
%         colsleft = [];
%         for i = candcols
%             Eltemp = eye(size(A));
%             nonbzeros = find(abs(B(:,i)) > thresh );
%             
%             if length(nonbzeros) > 1
%                 subteesubtorf = 0;
%                 tmps = transpose(combnk(nonbzeros,2));
%                 for rowcomb = tmps
%                     %                 rowcomb = nonbzeros(randperm(length(nonbzeros),2));
% 
%                     j1 = find(abs(Btemp(rowcomb(1),:)) > thresh);
%                     j2 = find(abs(Btemp(rowcomb(2),:)) > thresh);
%                     j1minusj2 = nonzeroacols(setdiff(j1,j2));
%                     j2minusj1 = nonzeroacols(setdiff(j2,j1));
%                     
%                     
%                     if isempty(intersect(j1minusj2, candcols(1:find(candcols==i,1)))) && isempty(intersect(j1minusj2, colsdone))
%                         subtractor = rowcomb(1);
%                         subtractee = rowcomb(2);
%                         if ~ismember(subtractor, subtractors)
%                             Eltemp = eye(size(A));
%                             Eltemp(subtractee,subtractor) = - B(subtractee,i)/B(subtractor,i);
%                             B = Eltemp * B;
%                             B(subtractee, i) = 0;
%                             A = Eltemp * A;
%                             Btemp = B(:,find(sum(abs(A))< thresh));
%                             fprintf(afile, '%s\n', strcat('A(', num2str(subtractee),',:) = A(',num2str(subtractee),',:) - A(', num2str(subtractor), ',:) * B(',num2str(subtractee),',',num2str(i), ')/B(',num2str(subtractor),',',num2str(i),');'));
%                             fprintf(afile, '%s\n', strcat('B(', num2str(subtractee),',:) = B(',num2str(subtractee),',:) - B(', num2str(subtractor), ',:) * B(',num2str(subtractee),',',num2str(i), ')/B(',num2str(subtractor),',',num2str(i),');'));
%                             subteesubtorf = 1;
%                         end
%                     elseif isempty(intersect(j2minusj1, candcols(1:find(candcols==i,1)))) && isempty(intersect(j2minusj1, colsdone))
%                         subtractor = rowcomb(2);
%                         subtractee = rowcomb(1);
%                         if ~ismember(subtractor, subtractors)
%                             Eltemp = eye(size(A));
%                             Eltemp(subtractee,subtractor) = - B(subtractee,i)/B(subtractor,i);
%                             B = Eltemp * B;
%                             B(subtractee, i) = 0;
%                             A = Eltemp * A;
%                             Btemp = B(:,find(sum(abs(A))< thresh));
%                             fprintf(afile, '%s\n', strcat('A(', num2str(subtractee),',:) = A(',num2str(subtractee),',:) - A(', num2str(subtractor), ',:) * B(',num2str(subtractee),',',num2str(i), ')/B(',num2str(subtractor),',',num2str(i),');'));
%                             fprintf(afile, '%s\n', strcat('B(', num2str(subtractee),',:) = B(',num2str(subtractee),',:) - B(', num2str(subtractor), ',:) * B(',num2str(subtractee),',',num2str(i), ')/B(',num2str(subtractor),',',num2str(i),');'));
%                             subteesubtorf = 1;
%                         end
%                     elseif length(j1) == length(j2)  && isempty(intersect(j1minusj2, colsdone))
%                         if j1 == j2
%                             subtractor = rowcomb(1);
%                             subtractee = rowcomb(2);
%                             if ~ismember(subtractor, subtractors)
%                                 Eltemp = eye(size(A));
%                                 Eltemp(subtractee,subtractor) = - B(subtractee,i)/B(subtractor,i);
%                                 B = Eltemp * B;
%                                 B(subtractee, i) = 0;
%                                 A = Eltemp * A;
%                                 Btemp = B(:,find(sum(abs(A))< thresh));
%                                 fprintf(afile, '%s\n', strcat('A(', num2str(subtractee),',:) = A(',num2str(subtractee),',:) - A(', num2str(subtractor), ',:) * B(',num2str(subtractee),',',num2str(i), ')/B(',num2str(subtractor),',',num2str(i),');'));
%                                 fprintf(afile, '%s\n', strcat('B(', num2str(subtractee),',:) = B(',num2str(subtractee),',:) - B(', num2str(subtractor), ',:) * B(',num2str(subtractee),',',num2str(i), ')/B(',num2str(subtractor),',',num2str(i),');'));
%                                 subteesubtorf = 1;
%                             end
%                         end
%                     end
%                     if subteesubtorf ==1
%                         break;
%                     end
%                 end
%             end
%             
%             nonbcolele = find(abs(B(:,i))> thresh);
%             if length(nonbcolele) ~= 1
%                 colsleft = [colsleft, i];
%             else
%                 subtractors=[subtractors, nonbcolele];
%                 colsdone = [colsdone, i];
%             end
%         end
%         if length(candcols) == length(colsleft)
%             if prod(candcols == colsleft)==1 && (retries <= 1)
%                 break;
%             end
%         end
%         retries = retries -1;
%         candcols = colsleft;
%     end
%     
%     rowperms = [];
%     colperms = [];
%     for i=colsdone
%         nonzerobcolrow = find(abs(B(:,i))>thresh);
%         if length(nonzerobcolrow)==1
%             rowperms=[rowperms, nonzerobcolrow];
%             colperms = [colperms,i];
%         end
%     end
%     if ~isempty(colperms) && ~isempty(rowperms) && length(unique(rowperms)) == length(rowperms)
%         A(:, colperms)=[];
%         B(:, colperms)=[];
%         A(rowperms,:)=[];
%         B(rowperms,:)=[];
%         perms = [colperms, rowperms];
%         extendedbasis(:,colperms) = [];
%         noofrows2rem = length(perms)/2;
%         fprintf(afile, '%s\n', strcat('col = ',mat2str(perms(1:noofrows2rem)),';'));
%         fprintf(afile, '%s\n', strcat('row = ', mat2str(perms(noofrows2rem +1:end)),';'));
%         fprintf(afile, '%s\n', 'A(:,col) = []; B(:,col) = []; A(row,:) = []; B(row,:) = [];');
%     end
%     %         turncolsinto0;    
%     fprintf("Removed  %d more z-columns and corresponding rows \n", size(colsdone,2));
%     fprintf("And then the matrix condition number is %f \n", cond(A));
%     
%     nonzeroacols = find(sum(abs(A)) < thresh);
%     candcols = nonzeroacols;
%     ext = ext - 1;
% end

fprintf(afile, '%s\n', 'end');
fclose(afile);
