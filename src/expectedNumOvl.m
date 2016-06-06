function ET = expectedNumOvl(m,S,Q)

SexpQtheta=sum(S,2);
c=1./SexpQtheta;
P=diag(c)*S*diag(expQtheta);
ET=m*P*Q;
