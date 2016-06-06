function grad = gradLogLikelihood(theta_hat,m,T,S,Q)

c=1./(S*exp(Q*theta_hat));
P=diag(c)*S*diag(exp(Q*theta_hat));
grad=T-m*P*Q;