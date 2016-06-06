function [rerun,param_ind_rerun] = testForHidden0ovl(pi,T,param_ind_init)

[~,N]=size(pi);
pow2=2.^(0:N-1);
piidx=pi*pow2';
T=T(piidx)';

% Find multiplicity of T
[m,T0]=hist(T,unique(T)); 
T0(m==0)=[];
m(m==0)=[];


idx=find(m>1);
nT=length(idx);
rerun=false;
param_ind_rerun=param_ind_init;
for i=1:nT
    pi_test=pi(T==T0(idx(i)),:);
    [M,~]=size(pi_test);
    for j=1:M-1
        for k=j+1:M
            if sum(pi_test(j,:).*pi_test(k,:))==sum(pi_test(j,:))
                % Remove pi_test(j,:) from model
                param_ind_rerun(pow2*pi_test(j,:)')=0;
                rerun=true;
            elseif sum(pi_test(j,:).*pi_test(k,:))==sum(pi_test(k,:))
                % Remove pi_test(k,:) from model
                param_ind_rerun(pow2*pi_test(k,:)')=0;
                rerun=true;
            end
        end
    end
end
    
