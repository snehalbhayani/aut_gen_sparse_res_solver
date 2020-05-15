thresh = 1e-6;
A2=A;
B2=B;
sizeofA2 = size(A2,1);


[sorted, soperms] =sort(sum(abs(A2) > thresh,1));
nullifiedcols = soperms(find(sorted<0));
soperms = soperms(find(sorted>0));
leftovercols = [];
% dividing the matrix into two parts such that one part is full rank but
% the other is not.
% temp = ceil(length(soperms)/2);
temp = sizeofA2;
fprintf("Using %d cols to eliminate", temp);
soperms = soperms(1:end);
rankofA2 = rank(A2);
s = soperms(1);
while length(soperms) > rankofA2    
    % Find cols with rows of non zero elements intersecting the given row
    % in non zero elements.
    candcols = [];
    for nonzerorow = transpose(find(abs(A2(:,s)) > thresh))
        candcols = transpose(union(candcols, transpose(find(abs(A2(nonzerorow,:)) > thresh))));
        candcols = candcols(find(candcols >= s));
    end
    noofnzins = length(candcols); 
    for k = 1:noofnzins
        broke = 0;
        candcols = setdiff(candcols, nullifiedcols);
        for comb = transpose(combnk(candcols,k))
            if rank(A2(:,comb)) == length(comb) - 1 && isempty(intersect(comb, nullifiedcols))
                nullvec = null(A2(:,comb));
%                 disp(nullvec);
%                 disp(" .... ................................. ");
                Eltemp = eye(sizeofA2);
                nnullvecr = find(abs(nullvec)>thresh,1);                
                fprintf(afile, '%s\n', strcat('nullvec=null(A(:,', mat2str(comb),'));'));
                
                fprintf(afile, '%s\n', strcat('B(:,', num2str(comb(nnullvecr)),') = B(:,',mat2str(comb),')*nullvec(:,1);'));
                for t = 1:length(comb)
                    if t ~= nnullvecr
                        Eltemp(comb(t),comb(nnullvecr)) = nullvec(t)/nullvec(nnullvecr);
                    end
                end
                A2 = A2 * Eltemp;
                B2 = B2 * Eltemp;
                nullifiedcols = [nullifiedcols, comb(nnullvecr)];
                
                [sorted, soperms] =sort(sum(abs(A2) > thresh,1));
                soperms = soperms(find(sorted>0));
                %             temp = ceil(length(soperms)/2);
                soperms = soperms(1:end);                
                s = soperms(find(soperms~=s,1));
                broke = 1;
                break;
            end
        end
    end
    
    if broke == 0 
        nextsind = find(soperms==s,1);
        if nextsind>=length(soperms)
            break;
        else
            s = soperms(nextsind+1);
        end
    end
    
end
fprintf(afile, '%s\n', strcat('A(:,', mat2str(nullifiedcols),') = 0;'));
A=A2;B=B2;
disp(sorted)
% for s=2:sizeofA1-rank(A2)
%     for comb = transpose(combnk(soperms,s))
%         if rank(A2(:,comb)) == s-1
%             nullvec = null(A2(:,comb));
%             Eltemp = eye(sizeofA1);
%             for t=2:s
%                 Eltemp(comb(t),comb(1)) = nullvec(t)/nullvec(1);
%             end
%             A2 = A2*Eltemp;
%             B2 = B2*Eltemp;
%             
%             [sorted, soperms] =sort(sum(abs(A2) > thresh,1));
%             soperms = soperms(find(sorted>0));
% %             temp = ceil(length(soperms)/2);
%             soperms = soperms(1:end);
%         end
%     end
% end