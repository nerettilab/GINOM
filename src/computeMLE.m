function [theta_hat,logL] = computeMLE(param_ind,m,T,S,Q)

N=log2(length(Q(:,1)));
nmod=sum(param_ind);

if nmod==2^N-1
    options=optimset('GradObj','on','Display','off');
    theta0=zeros(2^N-1,1);
    [theta_hat,logL]=fminunc(@(theta) neglogLikelihood(theta,m,T,S,Q),theta0,options);
    logL=-logL;
elseif nmod==0
    theta_hat=zeros(2^N-1,1);
    [logL,~]=logLikelihood(theta_hat,m,T,S,Q);
else
    options=optimset('GradObj','on','Display','off');
    A=[];
    b=[];
    Aeq=zeros(2^N-1-nmod,2^N-1);
    idx=find(param_ind==0);
    for i=1:length(idx)
        Aeq(i,idx(i))=1;
    end
    beq=zeros(2^N-1-nmod,1);
    lb=[];
    ub=[];
    nonlcon=[];
    theta0=zeros(2^N-1,1);
    [theta_hat,logL]=fmincon(@(theta) neglogLikelihood(theta,m,T,S,Q),theta0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    logL=-logL;
end


    