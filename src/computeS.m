function S = computeS(y,L,R,M)

n=length(y);
N=length(R);
S=zeros(n,2^N);
v=2.^(0:N-1);
for i=1:n
    for x=1:L-y(i)+1
        
        % Compute b
        b=computeb(x,y(i),R);
        
        % Convert b to S column index
        j=1+sum(v.*b);
        
        % Increment S
        S(i,j)=S(i,j)+f0(x,y(i),M);
        
    end
end
