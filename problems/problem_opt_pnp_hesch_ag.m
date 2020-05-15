function [ eqs, data0, eqs_data ] = problem_opt_pnp_hesch( data0 )

if nargin < 1 || isempty(data0)
    data0 = randi(50,35,1);
end


xx = create_vars(3);
vv = [1 xx(1) xx(2) xx(3) xx(1)^2 xx(1)*xx(2) xx(1)*xx(3) xx(2)^2 xx(2)*xx(3) xx(3)^2];
for iii = 1:3,
    for jjj = iii:3,
        for kkk = jjj:3,
            vv = [vv xx(iii)*xx(jjj)*xx(kkk)];
        end
    end
end
for iii = 1:3,
    for jjj = iii:3,
        for kkk = jjj:3,
            for lll = kkk:3
                vv = [vv xx(iii)*xx(jjj)*xx(kkk)*xx(lll)];
            end
        end
    end
end

f = vv*data0;
eqs = diff(f);
eqs = eqs(:);


if nargout == 3
    xx = create_vars(3+35);
    data = xx(4:end);
    vv = [1 xx(1) xx(2) xx(3) xx(1)^2 xx(1)*xx(2) xx(1)*xx(3) xx(2)^2 xx(2)*xx(3) xx(3)^2];
    for iii = 1:3,
        for jjj = iii:3,
            for kkk = jjj:3,
                vv = [vv xx(iii)*xx(jjj)*xx(kkk)];
            end
        end
    end
    for iii = 1:3,
        for jjj = iii:3,
            for kkk = jjj:3,
                for lll = kkk:3
                    vv = [vv xx(iii)*xx(jjj)*xx(kkk)*xx(lll)];
                end
            end
        end
    end

    f = vv*data(:);
    eqs_data = diff(f,1:3);
    eqs_data = eqs_data(:);

    
end

end