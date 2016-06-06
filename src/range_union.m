function out = range_union(first,second)

a=length(first);
firstl=zeros(1,a/2);
firstu=zeros(1,a/2);
for i=1:a/2
    firstl(i)=first(2*i-1);
    firstu(i)=first(2*i);
end

b=length(second);
secondl=zeros(1,b/2);
secondu=zeros(1,b/2);
for i=1:b/2
    secondl(i)=second(2*i-1);
    secondu(i)=second(2*i);
end

[l,u]=IntervalUnion([firstl secondl], [firstu secondu]);

n=length(u);
out=zeros(1,2*n);
for i=1:n
    out(2*i-1)=l(i);
    out(2*i)=u(i);
end

    