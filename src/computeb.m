function b = computeb(x,y,R)

N=length(R);
b=zeros(1,N);
for i=1:N
    Rext=[R{i}(:,1)-y+1 R{i}(:,2)];
    temp1=find(x>=Rext(:,1));
    temp2=find(x<=Rext(:,2));
    if ~isempty(temp1) && ~isempty(temp2)
        if temp2(1)==temp1(end) % If x lies within an interval of [R{i}(:,1)-y+1, R{i}(:,2)]
            b(i)=1;
        end
    end
end