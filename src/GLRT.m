function p = GLRT(param_ind_null,param_ind_full,m,T,S,Q)

nparam=length(param_ind_full);
diff=param_ind_full-param_ind_null;
if sum(diff<0)>0 || sum(diff==zeros(nparam,1))==nparam
    disp('Null parameter space must be a subset of the full parameter space');
    p='err';
else
    [~,logLnull] = computeMLE(param_ind_null,m,T,S,Q);
    [~,logLfull] = computeMLE(param_ind_full,m,T,S,Q);
    D=-2*(logLnull-logLfull);
    df=sum(diff);
    p=1-chi2cdf(D,df);
end
