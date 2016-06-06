function [Qrow,rownum] = getRowQ(b,Q)

N=length(b);
rownum=1+sum(2.^(0:N-1).*b);
Qrow=Q(rownum,:);