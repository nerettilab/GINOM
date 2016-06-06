function Ac = complement(Omega,A)

% A=[lower1, upper1, lower2, upper2,...]
% Omega is domain. In same format as A. Assume A \subset Omega

nOmega=length(Omega)/2;

Ac=[];
for i=1:nOmega
    
    x=A(A>=Omega(2*i-1) & A<=Omega(2*i));
    
    if isempty(x)
        
        Ac=[Ac Omega(2*i-1) Omega(2*i)];
        
    else
        
        flag1=x(1)==Omega(2*i-1);
        flag2=x(end)==Omega(2*i);

        if flag1 && ~flag2
            shift=[1 repmat([-1 1],1,length(x)/2-1)];
            Ac=[Ac x(2:end)+shift Omega(2*i)];
        elseif ~flag1 && flag2
            shift=[repmat([-1 1],1,length(x)/2-1) -1];
            Ac=[Ac Omega(2*i-1) x(1:end-1)+shift];
        elseif flag1 && flag2
%             shift=[1 repmat([-1 1],1,length(x)/2-2) -1];
            if length(x)~=2
                shift=repmat([1 -1],1,length(x)/2-1);
                Ac=[Ac x(2:end-1)+shift];
            end
        else 
            shift=repmat([-1 1],1,length(x)/2);
            Ac=[Ac Omega(2*i-1) x+shift Omega(2*i)];
        end
        
    end
    
end