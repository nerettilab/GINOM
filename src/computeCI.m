function CI = computeCI(theta_hat,param_ind,m,S,Q)

% Compute Fisher info
c=1./(S*exp(Q*theta_hat));
P=diag(c)*S*diag(exp(Q*theta_hat));
I=Q'*diag(m*P)*Q-(P*Q)'*diag(m)*(P*Q);

% Compute 95% confidence intervals
Iinv=inv(I);
CI(:,1)=theta_hat-1.96*sqrt(diag(Iinv));
CI(:,2)=theta_hat+1.96*sqrt(diag(Iinv));
CI=CI(param_ind==1,:);
    