function [theta_hat,param_ind_model] = stepwise(penalty,param_ind,m,T,S,Q)

K=length(param_ind);
nmax=sum(param_ind);
param_ind_model=zeros(K,1);
param_ind_remaining=param_ind;
[IC,~,~] = computeIC(penalty,param_ind_model,m,T,S,Q);

while 1<2 % Sometimes "break" may never be reached. Careful of infinite loop.
    
    nmod=sum(param_ind_model);
    nrem=sum(param_ind_remaining);
    
    if nmod==0 % If there are no model terms currently 
        
        % Try to add a model term
        IC_add=zeros(nrem,1);
        idx=find(param_ind_remaining==1);
        for i=1:nrem
            param_ind_try=param_ind_model;
            param_ind_try(idx(i))=1;
            [IC_add(i),~,~]=computeIC(penalty,param_ind_try,m,T,S,Q);
        end
        if sum(IC_add<IC)==0 % If the addition of no model term lowers the IC, the null model is best model, and break
            theta_hat=[];
%             logL=logLnull;
            break;
        else % If the addition of at least one model term lowers the IC, add the one that yields the lowest value 
            [IC,idx_add]=min(IC_add);
            param_ind_model(idx(idx_add))=1;
            param_ind_remaining(idx(idx_add))=0;
        end
        
    else % If there is at least one model term currently
        
        % Try to delete a model term
        IC_del=zeros(nmod,1);
        idx=find(param_ind_model==1);
        for i=1:nmod
            param_ind_try=param_ind_model;
            param_ind_try(idx(i))=0;
            [IC_del(i),~,~]=computeIC(penalty,param_ind_try,m,T,S,Q);
        end
        
        if sum(IC_del<IC)==0 % If the deletion of no model term lowers the IC, 
            
            if nmod<nmax % If nmod < nmax, try to add a model term
                IC_add=zeros(nrem,1);
                idx=find(param_ind_remaining==1);
                for i=1:nrem
                    param_ind_try=param_ind_model;
                    param_ind_try(idx(i))=1;
                    [IC_add(i),~,~]=computeIC(penalty,param_ind_try,m,T,S,Q);
                end
                if sum(IC_add<IC)==0 % If the addition of no model term lowers the IC, the current model is best model, and break
                    [~,~,theta_hat]=computeIC(penalty,param_ind_model,m,T,S,Q);
                    break;
                else % If the addition of at least one model term lowers the IC, add the one that yields the lowest value
                    [IC,idx_add]=min(IC_add);
                    param_ind_model(idx(idx_add))=1;
                    param_ind_remaining(idx(idx_add))=0;
                end
            else % If nmod==nmax, the full model is the best model, and break
                [~,~,theta_hat]=computeIC(penalty,param_ind_model,m,T,S,Q);
                break;
            end
                
        else % If the deletion of at least one model term lowers the IC, remove the one that the removal of which would yield the lowest IC
            
            [IC,idx_rmv]=min(IC_del);
            param_ind_model(idx(idx_rmv))=0;
            param_ind_remaining(idx(idx_rmv))=1;

        end
        
    end
    
end
