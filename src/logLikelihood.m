function [ell,gradell] = logLikelihood(theta,m,T,S,Q)

Qtheta=Q*theta;
expQtheta=exp(Qtheta);
SexpQtheta=S*expQtheta;
ell=T*theta-m*log(SexpQtheta);

c=1./SexpQtheta;
P=diag(c)*S*diag(expQtheta);
gradell=T-m*P*Q;

% c=S*exp(Q*theta);
% c=1./c;
% cvals=c(yind);
% % logfxy=sum(log(cvals))+sum(log(f0vals))+sum(rq*theta);
% logfxy=sum(log(cvals))+sum(rq*theta);