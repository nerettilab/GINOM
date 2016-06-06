function [Qcol,colnum] = getColQ(pi,Q)

N=length(pi);
colnum=sum(2.^(0:N-1).*pi);
Qcol=Q(:,colnum);