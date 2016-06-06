function [theta_hat,param_ind_model] = modelSelectFromAll(penalty,param_ind,m,T,S,Q)

% keyboard;
K=length(param_ind);
nmax=sum(param_ind);
combos=dec2bin(0:2^nmax-1);
[a,b]=size(combos);
idx=zeros(a,b);
for i=1:a      
    for j=1:b
        idx(i,j)=str2double(combos(i,j));
    end
end

% Compute IC for all model combinations
param_ind_test=zeros(a,K);
IC=zeros(a,1);
theta_hat_test=zeros(K,a);
h=waitbar(0,'Fitting all possible models');
for i=1:a
    param_ind_test(i,param_ind==1)=idx(i,:);
    [IC(i),~,theta_hat_test(:,i)] = computeIC(penalty,param_ind_test(i,:),m,T,S,Q);
    waitbar(i/a);
end
close(h);

% Pick model with minimum IC
[~,minidx]=min(IC);
theta_hat=theta_hat_test(:,minidx);
param_ind_model=param_ind_test(minidx,:)';