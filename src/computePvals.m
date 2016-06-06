function pvals = computePvals(param_ind,m,T,S,Q)

model_idx=find(param_ind==1);
nmodel=length(model_idx);
pvals=zeros(nmodel,1);
for i=1:nmodel
    param_ind_null=param_ind;
    param_ind_null(model_idx(i))=0;
    pvals(i) = GLRT(param_ind_null,param_ind,m,T,S,Q);
end