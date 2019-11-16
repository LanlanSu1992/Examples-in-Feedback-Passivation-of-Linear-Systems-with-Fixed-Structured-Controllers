function [J,C] = modified_jury(coeff)


J = [coeff;flipdim(coeff,2)];
typ = class(coeff);
n = length(coeff)-1;

if strcmp(typ,'sym')
    for i=3:2:(2*n+1)
        try
            alph1 = J(i-1,1);
            alph2=J(i-2,1);
        catch
            disp('Your polynomial seems to be critical')
            rethrow(lasterror);
            break;
        end
        newrow_1 = J(i-2,:)*alph2-alph1*J(i-1,:);
        %newrow = simplify(newrow_1);
        newrow=newrow_1;
        J = [J ; newrow ; 
                [flipdim(newrow(1:end-(i-1)/2),2) , zeros(1,(i-1)/2)]
                                    ];
    end
else
    for i=3:2:(2*n+1)
        try
            alph1 = J(i-1,1);
            alph2=J(i-2,1);
        catch
            disp('Your polynomial seems to be critical')
            rethrow(lasterror);
            break;
        end
        newrow = J(i-2,:)*alph2-alph1*J(i-1,:);
        J = [J ; newrow ; 
                [flipdim(newrow(1:end-(i-1)/2),2) , zeros(1,(i-1)/2)]
                                    ];
    end
end

J = J(1:end-1,:);
C = J(1:2:end,1);


