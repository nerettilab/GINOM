function preprocessData(filename)

disp('Reading input text file');
tic
timeidx=0;

fid=fopen(filename);
C=textscan(fid,'%s','Delimiter','|');

reference=C{1}{1};
query=C{1}{2};
mappability=C{1}{3};
chr_lengths=C{1}{4};
data_path=C{1}{5};
output_path=C{1}{6};
output_name=C{1}{7};
chr=textscan(C{1}{8},'%d','Delimiter',',');
chr=chr{1};

% Keep track of computation time
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

disp('Loading data');
tic

% Add data path to working directory
addpath(data_path);

% Load reference interval sets, R 
load(reference); 

% Number of reference interval sets
N=length(R); %#ok<NODEF>

% Load chromosome lengths, L
load(chr_lengths);

% keyboard;

% Concatenate ref interval start/end indices (put one chromosome after 
% another)
Rorig=R;
clear R;
R=struct([]);
for i=1:N
    R(i).Chr=[];
    R(i).Start=[];
    R(i).End=[];
    l=0;
    for j=1:length(chr)
        chr_idx=strcmp(L.Chr{chr(j)},Rorig(i).Chr);
        R(i).Chr=[R(i).Chr; Rorig(i).Chr(chr_idx)];
        R(i).Start=[R(i).Start; Rorig(i).Start(chr_idx)+l];
        R(i).End=[R(i).End; Rorig(i).End(chr_idx)+l];
        l=l+L.Len(chr(j));
    end
end
clear Rorig;

% Load mappability data, 'Map'.
load(mappability);

% Concatenate mappability start/end indices (put one chromosome after 
% another)
MapOrig=Map; %#ok<NODEF>
clear Map;
l=0;
Map.Chr=[];
Map.Start=[];
Map.End=[];
for i=1:length(chr)
    chr_idx=strcmp(L.Chr{chr(i)},MapOrig.Chr);
    Map.Chr=[Map.Chr; MapOrig.Chr(chr_idx)];
    Map.Start=[Map.Start; MapOrig.Start(chr_idx)+l];
    Map.End=[Map.End; MapOrig.End(chr_idx)+l];
    l=l+L.Len(chr(i));
end
clear MapOrig;

% Load query intervals, Q
load(query);

% Concatenate query data, Q, as above (and sort on each chr). 
Qorig=Q; %#ok<NODEF>
clear Q;
l=0;
Q.Chr=[];
Q.Start=[];
Q.End=[];
for i=1:length(chr)
    
    chr_idx=strcmp(L.Chr{chr(i)},Qorig.Chr);
    
    Q.Chr=[Q.Chr; Qorig.Chr(chr_idx)];
    
    [temp1,idx]=sort(Qorig.Start(chr_idx)+l); 
    Q.Start=[Q.Start; temp1]; 
    
    temp2=Qorig.End(chr_idx)+l;
    temp2=temp2(idx);
    Q.End=[Q.End; temp2];
    
    l=l+L.Len(chr(i));
    
end
clear Qorig;

% Assign qx to be query interval start indices and qy to be the query 
% interval lengths
qx=Q.Start;
qy=Q.End-Q.Start+1;
clear Q;

clear chr_idx temp1 temp2;

% Keep track of computation time
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

% keyboard;

% Here, all reference interval sets are sorted and unioned.
disp('Sorting, unioning R');
tic
[Map,R]=preprocessR(Map,R);
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

% keyboard;

% Create M structure
disp('Creating M structure');
tic
M=formMstruct(R,Map);
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

% Compute unique values of qy and their multiplicities m
disp('Computing unique values of y');
tic
[m,y]=hist(qy,unique(qy)); 
y(m==0)=[];
m(m==0)=[]; %#ok<NASGU>
% m=m'; 
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

% Compute S
disp('Computing S');
tic
S=computeS_f0uniform(M,y); %#ok<NASGU>
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

% Compute Q
disp('Computing Q');
tic
Q=createQ(N);
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

% Compute r for each query interval.
disp('Computing r for each query interval')
tic
n=length(qx);
rq=zeros(n,2^N-1);
for i=1:n
    rq(i,:)=getRowQ(computeb(qx(i),qy(i),R),Q);
end
T=sum(rq); %#ok<NASGU>
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

% theta index scheme
disp('Computing set indexing scheme pi')
tic
pow2=2.^(0:N-1);
combos=dec2bin(1:2^N-1);
pi=zeros(2^N-1,N);
npi=zeros(2^N-1,1);
for i=1:2^N-1        % (there's probably an efficient way to do this)
    for j=1:N
        pi(i,j)=str2double(combos(i,j));
    end
    npi(i)=pow2*pi(i,:)';
end
pi=pi(npi,:);
b=[zeros(1,N);pi]; %#ok<NASGU> % Binary indexing scheme
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

% keyboard;

disp('Saving output variables')
tic
save([output_path '\' output_name],'S','Q','T','m','pi','b');
timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');
disp(['Total elapsed time is ' num2str(sum(time)) ' seconds.']) 

