function [IC,logL,theta_hat] = computeIC(penalty,param_ind,m,T,S,Q)

M=sum(param_ind);
[theta_hat,logL]=computeMLE(param_ind,m,T,S,Q);
if strcmp(penalty,'AIC')
    IC=2*M-2*logL;
elseif strcmp(penalty,'SBC')
    IC=M*log(sum(m))-2*logL;
end
