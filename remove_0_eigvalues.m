%% Stage 1 Removal of 0 eigen values through same row and column
    
    zeigvalindx = [];
    temp = A(end - sizeOfReducedCs +1:end,:);    
    zcol = find(sum(abs(A(:,1:end-sizeOfReducedCs))) == 0,1);
    while(~isempty(zcol))
        A(:,zcol) = [];
        A(zcol,:) = [];
        B(zcol,:) = [];
        B(:,zcol) = [];
        temp(:,zcol) = [];
        extendedbasis(:,zcol) = [];        
        zeigvalindx = [zeigvalindx, zcol];        
        zcolnext = find(sum(abs(A(:,1:end-sizeOfReducedCs))) == 0);
        zcol = zcolnext(find(zcolnext>=zcol,1));

        while(~isempty(zcol))
            A(:,zcol) = [];
            A(zcol,:) = [];
            B(zcol,:) = [];
            B(:,zcol) = [];
            temp(:,zcol) = [];
            extendedbasis(:,zcol) = [];
            zeigvalindx = [zeigvalindx, zcol];
            zcolnext = find(sum(abs(A(:,1:end-sizeOfReducedCs))) == 0);
            zcol = zcolnext(find(zcolnext>=zcol,1));            
        end
        zcol = find(sum(abs(A(:,1:end-sizeOfReducedCs))) == 0,1);        
    end
    

%% Stage 2: Trying to remove the rows-columns that can be removed so that we can remove
%     % more 0 eigen values
%     % 1. Either we scan the whole A matrix for columns 
%     % 2. Or we skip the last 'sizeOfReducedCs' columns of matrix.
    colpermutation=[];
    rowpermutation=[];
    fullscanning = 1;
    temprange = size(A,1)-sizeOfReducedCs+1;
    temp = A(temprange:end,:);
    if fullscanning == 0
        zcolsina = find(sum(abs(temp(:,1:end-sizeOfReducedCs))) == 0);
    else
        zcolsina = find(sum(abs(A(:,1:end))) == 0);
    end
    it = 2;
    zcola = zcolsina(find(sum(abs(B(:,zcolsina)) == sum(abs(B(:,zcolsina)))) == 1,1));
%     zcola = find(zcolsina,1);
    zcolb = find(B(:,zcola) ~= 0, 1);
    
    while(~isempty(zcola))
        A(:,zcola) = [];
        A(zcolb,:) = [];
        B(zcolb,:) = [];
        B(:,zcola) = [];
        colpermutation = [colpermutation, zcola];
        rowpermutation = [rowpermutation, zcolb];
        temp = A(temprange:end,:);
        extendedbasis(:,zcola) = [];
        if fullscanning == 0
            zcolsina = find(sum(abs(A(:,1:end-sizeOfReducedCs))) == 0);
        else
            zcolsina = find(sum(abs(A(:,1:end))) == 0);
        end
        zcolnext = find(sum(abs(B(:,zcolsina)) == sum(abs(B(:,zcolsina)))) == 1);
        zcola = zcolsina(zcolnext(find(zcolsina(zcolnext)>=zcola, 1)));
%         zcola = find(zcolsina,1);
        zcolb = find(B(:,zcola) ~= 0, 1);
        while(~isempty(zcola))
            A(:,zcola) = [];
            A(zcolb,:) = [];
            B(zcolb,:) = [];
            B(:,zcola) = [];
            colpermutation = [colpermutation, zcola];
            rowpermutation = [rowpermutation, zcolb];
            temp = A(temprange:end,:);
            extendedbasis(:,zcola) = [];
            
            if fullscanning == 0
                zcolsina = find(sum(abs(A(:,1:end-sizeOfReducedCs))) == 0);
            else
                zcolsina = find(sum(abs(A(:,1:end))) == 0);
            end
            zcolnext = find(sum(abs(B(:,zcolsina)) == sum(abs(B(:,zcolsina)))) == 1);
            zcola = zcolsina(zcolnext(find(zcolsina(zcolnext)>=zcola, 1)));
%             zcola = find(zcolsina,1);
            zcolb = find(B(:,zcola) ~= 0, 1);
        end
        
        if fullscanning == 0
            zcolsina = find(sum(abs(temp(:,1:end-sizeOfReducedCs))) == 0);
        else
            zcolsina = find(sum(abs(A(:,1:end))) == 0);
        end
        zcola = zcolsina(find(sum(abs(B(:,zcolsina)) == sum(abs(B(:,zcolsina)))) == 1,1));
%         zcola = find(zcolsina,1);
        zcolb = find(B(:,zcola) ~= 0, 1);
        
    end