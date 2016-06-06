function [M,Rnew] = preprocessR(Map,R)

N=length(R);

% Mappability
nM=length(Map.Start);
M(:,1)=Map.Start;
M(:,2)=Map.End;

% Reorganize for future set operations
Omega=zeros(1,2*nM);
Omega(1:2:2*nM-1)=M(:,1);
Omega(2:2:2*nM)=M(:,2);

% Preprocess regions of interest
Rnew{N}=[];
nR=zeros(N,1);
r{N}=[];
for i=1:N
    
    % Sort them
    [Rnew{i}(:,1),idx]=sort(R(i).Start);
    Rnew{i}(:,2)=R(i).End(idx);
    
    % Compute union of the intervals so that the interval limits are 
    % disjoint
    [temp(1,:),temp(2,:)]=IntervalUnion(Rnew{i}(:,1)',Rnew{i}(:,2)');
    Rnew{i}=temp';
    clear temp;
    
    % Merge adjacent intervals
    Rnew{i}=MergeAdjacentIntervals(Rnew{i});
    [nR(i),~]=size(Rnew{i});
    
    % Intersect R{i} with mappable region
    r{i}=zeros(1,2*nR(i));
    r{i}(1:2:2*nR(i)-1)=Rnew{i}(:,1);
    r{i}(2:2:2*nR(i))=Rnew{i}(:,2);
    r{i}=range_intersection(Omega,r{i});
    
    % Redefine Rnew after the above intersection
    nR(i)=length(r{i})/2;
    Rnew{i}=zeros(nR(i),2);
    Rnew{i}(:,1)=r{i}(1:2:2*nR(i)-1);
    Rnew{i}(:,2)=r{i}(2:2:2*nR(i));
    
    % Merge adjacent intervals
    Rnew{i}=MergeAdjacentIntervals(Rnew{i});
    [nR(i),~]=size(Rnew{i});

end