function intervals = MergeAdjacentIntervals(G)

% G is a Nx2 matrix representing the starting and ending point of each
% iterval.

% Merge adjacent intervals
d=G(2:end,1)-G(1:end-1,2);
neq=find(d~=1);
n=length(neq);
intervals=zeros(n+1,2);
intervals(1,1)=G(1,1);
for i=2:n+1
    intervals(i-1,2)=G(neq(i-1),2);
    intervals(i,1)=G(neq(i-1)+1,1);
end
intervals(n+1,2)=G(end,2);