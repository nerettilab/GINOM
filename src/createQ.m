function Q = createQ(N)

Q=[0;1];
A=[0 1];
if N>1
    for i=1:N-1
        first=[Q;Q];
        [r,c]=size(first);
        newcol=reshape(repmat(A(:).',2^i,1),1,[])';
        Q=[first,newcol,reshape(first(:).*repmat(newcol,c,1),r,c)];
    end
end