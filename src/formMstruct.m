function M = formMstruct(R,Map)

% Map is a list of disjoint mappable intervals, nMap x 2
% R is a cell array of reference intervals, each nR(i) x 2. Each R{i} is
% sorted and consists of disjoint intervals contained entirely within Map. 

N=length(R);

% Reorganize Map for future set operations
nM=length(Map(:,1));
map=zeros(1,2*nM);
map(1:2:2*nM-1)=Map(:,1)';
map(2:2:2*nM)=Map(:,2)';

% keyboard;

% Reorganze R for future set operations
r{N}=[];
nR=zeros(N,1);
for i=1:N
    nR(i)=length(R{i}(:,1));
    r{i}=zeros(1,2*nR(i));    
    r{i}(1:2:2*nR(i)-1)=R{i}(:,1)';
    r{i}(2:2:2*nR(i))=R{i}(:,2)';
end

% Handle binary indexing 
pow2=2.^(0:N-1);
combos=dec2bin(0:2^N-1)';
b=zeros(N,2^N);
nb=zeros(1,2^N);
for j=1:2^N        % (there's probably an efficient way to do this)
    for i=1:N
        b(i,j)=str2double(combos(i,j));
    end
    nb(j)=1+pow2*b(:,j);
end
b=b(:,nb);

% Create disjointized R, say D, where D{1},D{2},...,D{2^N} is indexed by 
% binary index scheme b. 

D{2^N}=[];
nD=zeros(2^N,1);
for i=1:2^N
    
    if i==1 % If b(:,i)==[0;0;0]
    
        % Union all r and form d{1} as the complement map \ union r
        union=r{1};
        if N>1
            for j=2:N
                union=range_union(union,r{j});
            end
        end
        d=complement(map,union);
        
    elseif i==2^N % If b(:,i)==[1;1;1]
        
        % Form d{2^N} as the intersection of all r
        intersect=r{1};
        if N>1
            for j=2:N
                intersect=range_intersection(intersect,r{j});
            end
        end
        d=intersect;
        
    else
        
        union_idx=find(b(:,i)==0);
        intersect_idx=find(b(:,i)==1);
        
        intersect=r{intersect_idx(1)};
        if length(intersect_idx)>1
            for j=2:length(intersect_idx)
                intersect=range_intersection(intersect,r{intersect_idx(j)});
            end
        end

        union=r{union_idx(1)};
        if length(union_idx)>1
            for j=2:length(union_idx)
                union=range_union(union,r{union_idx(j)});
            end
        end

        d=range_intersection(intersect,complement(map,union));
        
    end
    
    % Form D{i}
    nD(i)=length(d)/2;
    D{i}(:,1)=d(1:2:2*nD(i)-1)';
    D{i}(:,2)=d(2:2:2*nD(i))';
        
end

% Make Dtot as one big list of all D intervals, indexed by D_nb
% Dtot=zeros(sum(nD),2);
% D_nb=zeros(sum(nD),1);
Dtot=[];
D_nb=[];
for i=1:2^N
    Dtot=[Dtot;D{i}];
    D_nb=[D_nb;i*ones(nD(i),1)];
end

% Sort Dtot (and D_nb accordingly)
[Dtot_sort,idx_sort]=sort(Dtot(:,1));
Dtot_sort(:,2)=Dtot(idx_sort,2);
D_nb_sort=D_nb(idx_sort);

% Form M structure
M=struct();
for i=1:nM
    
    imin=find(Dtot_sort(:,1)==Map(i,1));
    imax=find(Dtot_sort(:,2)==Map(i,2));
    
    k=0;
    M(i).b=zeros(N,imax-imin+1);
    M(i).w=zeros(1,imax-imin+1);
    for j=imin:imax
        k=k+1;
        M(i).b(:,k)=b(:,D_nb_sort(j));
        M(i).w(k)=Dtot_sort(j,2)-Dtot_sort(j,1)+1;
    end
    
end


% keyboard;


























% Convert b to binary representation to form matrix M.b.
