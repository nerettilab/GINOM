function S = computeS_f0uniform(M,y)

N=length(M(1).b(:,1)); % Number of reference sets
S=zeros(length(y),2^N);
for i=1:length(M)
    S=incrementS(S,y,M(i).b,M(i).w);
end

end

function S = incrementS(S,y,b,w)

imax=length(w);
ymax=max(y);
pow2=2.^(0:(length(b(:,1))-1));
for i=1:imax
    nbi=1+pow2*b(:,i);
    S=updateA(w(i),S,y,nbi);
    if i==imax
        break % return S
    end
    for j=(i+1):imax
        if j==i+1
            L=0;
        else
            L=sum(w((i+1):(j-1)));
        end
        if L+2>ymax
            break
        end
%         keyboard;
        bmax=max(b(:,i:j),[],2); % check this line
        nbmax=1+pow2*bmax;
        S=updateB(w(i),w(j),L,S,y,nbmax);
    end
end

end

function S = updateA(w,S,y,nb)

for i=1:length(y)
    if y(i)>w
        break
    end
    S(i,nb)=S(i,nb)+w-y(i)+1;
end

end

function S = updateB(w1,w2,L,S,y,nb)

low=L+2;
hi=L+w1+w2;
min12=min(w1,w2);
for i=1:length(y)
    yy=y(i);
    if yy<low
        continue
    elseif yy>hi
        break
    else
        S(i,nb)=S(i,nb)+min([yy-low+1,hi-yy+1,min12]);
    end
end

end
        