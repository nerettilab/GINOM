function [ell,gradell] = neglogLikelihood(theta,m,T,S,Q)

Qtheta=Q*theta;
expQtheta=exp(Qtheta);
SexpQtheta=S*expQtheta;
ell=T*theta-m*log(SexpQtheta);
ell=-ell;

c=1./SexpQtheta;
P=diag(c)*S*diag(expQtheta);
gradell=T-m*P*Q;
gradell=-gradell;