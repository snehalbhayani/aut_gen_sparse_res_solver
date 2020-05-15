function [ eqs, data0, eqs_data ] = problem_gen_26_v2( data0 )

xx = create_vars(3);


if nargin < 1 || isempty(data0)
    ll1 = randi(20,5,6);
    ll1 = [ll1;-(ll1(1,:).*ll1(4,:)+ll1(2,:).*ll1(5,:))];
    ll1(1:5,:)=ll1(1:5,:).*repmat(ll1(3,:),5,1);
    ll2 = randi(20,5,6);
    ll2 = [ll2;-(ll2(1,:).*ll2(4,:)+ll2(2,:).*ll2(5,:))];
    ll2(1:5,:)=ll2(1:5,:).*repmat(ll2(3,:),5,1);
    pp = 30097;
    ll1n2 = sum(ll1(1:3,:).^2);
    ll2n2 = sum(ll2(1:3,:).^2);
    R = 2*(xx*xx'-[0 -xx(3) xx(2);xx(3) 0 -xx(1);-xx(2) xx(1) 0])+(1-xx'*xx)*eye(3);
    M = zeros(30,84);
    %eq0 = [];
    Fk(5,3)=multipol;
    iddes = nchoosek(1:5,3);
    count2 = 0;
    for kkk = 1:3,
        x1 = cross(ll1(4:6,kkk),ll1(1:3,kkk));
        d1 = ll1(1:3,kkk);
        n1 = ll1n2(kkk);
        x2 = cross(ll2(4:6,kkk),ll2(1:3,kkk));
        d2 = ll2(1:3,kkk);
        n2 = ll2n2(kkk);
        
        id1 = setdiff(1:6,kkk);
        count = 0;
        for iii = id1,
            count = count + 1;
            q1 = ll1(1:3,iii);
            q1p = ll1(4:6,iii);
            q2 = ll2(1:3,iii);
            q2p = ll2(4:6,iii);
            
            tmp = n1*n2*q2'*R*q1p+q2'*(n2*R*[0 -x1(3) x1(2);x1(3) 0 -x1(1);-x1(2) x1(1) 0]-n1*[0 -x2(3) x2(2);x2(3) 0 -x2(1);-x2(2) x2(1) 0]*R)*q1+n1*n2*q1'*R'*q2p;
            cc = coeffs(tmp);
            mm = monomials(tmp);
            Fk(count,1)=multipol(mod(cc,pp),mm);
            %Fk(count,1)=multipol(cc,mm);
            tmp = n2*q2'*R*[0 -d1(3) d1(2);d1(3) 0 -d1(1);-d1(2) d1(1) 0]*q1;
            cc = coeffs(tmp);
            mm = monomials(tmp);
            Fk(count,2)=multipol(mod(cc,pp),mm);
            %Fk(count,2)=multipol(cc,mm);
       
            tmp = -n1*q2'*[0 -d2(3) d2(2);d2(3) 0 -d2(1);-d2(2) d2(1) 0]*R*q1;
            cc = coeffs(tmp);
            mm = monomials(tmp);
            Fk(count,3)=multipol(mod(cc,pp),mm);
            %Fk(count,3)=multipol(cc,mm);
        end
        
        
        for jjj = 1:size(iddes,1),
            eq = det(Fk(iddes(jjj,:),:));
            cc = mod(round(coeffs(eq)),pp);
            %mm = monomials(eq);
            %cc = coeffs(eq);
            count2 = count2 + 1;
            %eq0 = [eq0;multipol(cc,mm)];
            %disp(length(cc));
            M(count2,:)=cc;
        end
        
        
        
    end
    %M = polynomials2matrix(eq0);
    M = M([1:10 17 18 19 20 30],:);
    
    data0 = M(:);
end

e1 = [6 5 5 5 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
e2 = [0 1 0 0 2 1 1 0 0 0 3 2 2 1 1 1 0 0 0 0 4 3 3 2 2 2 1 1 1 1 0 0 0 0 0 5 4 4 3 3 3 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 6 5 5 4 4 4 3 3 3 3 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 ];
e3 = [0 0 1 0 0 1 0 2 1 0 0 1 0 2 1 0 3 2 1 0 0 1 0 2 1 0 3 2 1 0 4 3 2 1 0 0 1 0 2 1 0 3 2 1 0 4 3 2 1 0 5 4 3 2 1 0 0 1 0 2 1 0 3 2 1 0 4 3 2 1 0 5 4 3 2 1 0 6 5 4 3 2 1 0 ];
M = reshape(data0,15,84);
vv = xx(1).^e1.*xx(2).^e2.*xx(3).^e3;
eqs = M*vv(:);

if nargout == 3
    xx = create_vars(3+15*84);
    %     data = xx(4:end);
    %
    ll1 = reshape(xx(4:33),5,6);
    ll1 = [ll1;-(ll1(1,:).*ll1(4,:)+ll1(2,:).*ll1(5,:))];
    ll1(1:5,:)=ll1(1:5,:).*repmat(ll1(3,:),5,1);
    ll2 = reshape(xx(34:63),5,6);
    ll2 = [ll2;-(ll2(1,:).*ll2(4,:)+ll2(2,:).*ll2(5,:))];
    ll2(1:5,:)=ll2(1:5,:).*repmat(ll2(3,:),5,1);
    ll1n2 = sum(ll1(1:3,:).^2);
    ll2n2 = sum(ll2(1:3,:).^2);
    xx = xx(1:3);
    R = 2*(xx*xx'-[0 -xx(3) xx(2);xx(3) 0 -xx(1);-xx(2) xx(1) 0])+(1-xx'*xx)*eye(3);
    M = zeros(30,84);
    %eq0 = [];
    Fk(5,3)=multipol;
    iddes = nchoosek(1:5,3);
    count2 = 0;
    for kkk = 1:3,
        x1 = cross(ll1(4:6,kkk),ll1(1:3,kkk));
        d1 = ll1(1:3,kkk);
        n1 = ll1n2(kkk);
        x2 = cross(ll2(4:6,kkk),ll2(1:3,kkk));
        d2 = ll2(1:3,kkk);
        n2 = ll2n2(kkk);
        
        id1 = setdiff(1:6,kkk);
        count = 0;
        for iii = id1,
            count = count + 1;
            q1 = ll1(1:3,iii);
            q1p = ll1(4:6,iii);
            q2 = ll2(1:3,iii);
            q2p = ll2(4:6,iii);
            

            tmp = n1*n2*q2'*R*q1p+q2'*(n2*R*[0 -x1(3) x1(2);x1(3) 0 -x1(1);-x1(2) x1(1) 0]-n1*[0 -x2(3) x2(2);x2(3) 0 -x2(1);-x2(2) x2(1) 0]*R)*q1+n1*n2*q1'*R'*q2p;
            cc = coeffs(tmp);
            mm = monomials(tmp);
            Fk(count,1)=multipol(cc, mm);
            %Fk(count,1)=multipol(cc,mm);
            tmp = n2*q2'*R*[0 -d1(3) d1(2);d1(3) 0 -d1(1);-d1(2) d1(1) 0]*q1;
            cc = coeffs(tmp);
            mm = monomials(tmp);
            Fk(count,2)=multipol(cc,mm);
            %Fk(count,2)=multipol(cc,mm);
       
            tmp = -n1*q2'*[0 -d2(3) d2(2);d2(3) 0 -d2(1);-d2(2) d2(1) 0]*R*q1;
            cc = coeffs(tmp);
            mm = monomials(tmp);
            Fk(count,3)=multipol(cc,mm);
            %Fk(count,3)=multipol(cc,mm);
        end
        
        
        for jjj = 1:size(iddes,1),
            eq = det(Fk(iddes(jjj,:),:));
            cc = coeffs(eq);
            %mm = monomials(eq);
            %cc = coeffs(eq);
            count2 = count2 + 1;
            %eq0 = [eq0;multipol(cc,mm)];
            %disp(length(cc));
            M(count2,:)=cc;
        end
        
        
        
    end
    %M = polynomials2matrix(eq0);
    M = M([1:10 17 18 19 20 30],:);
    
    data = M(:);
    
    
    M = reshape(data,15,84);
    vv = xx(1).^e1.*xx(2).^e2.*xx(3).^e3;
    eqs_data = M*vv(:);
    
    totalsyms = 3 + 15*84;
    noofvars = 3;
    ipparams = strcat('a',num2str(1));
    for i = 2:totalsyms
        if i > noofvars
            ipparams = strjoin({ipparams, ',c', num2str(i-noofvars)}, '');
        else
            ipparams = strjoin({ipparams, ',a', num2str(i)}, '');
        end
    end
    
    fileID = fopen('Eqs_problem_gen_26_v2.m','w');
    fprintf(fileID, '%s', strjoin({'function eqs = retrieve_eqs(', ipparams, ') '},''));
    fprintf(fileID, '\n');
    for i = 1:size(eqs,1)
        fprintf(fileID, '%s', strjoin({'eqs(', num2str(i),') = ',char(eqs(i,1), true, [], true, noofvars), ';'},''));
        fprintf(fileID, '\n');
    end
    fprintf(fileID, '\n');
    fprintf(fileID, '%s', 'end');
    fclose(fileID);
end
