function [card1,card2,card12] = cardinality_intersection(R1,R2,chr)

% chr=19; % Try 18 too
% chr=1:24;

addpath('C:\Users\Darshan\Documents\MATLAB\Transposable Elements\TE Insertion Bias\Data');
load Chr_lengths_human

R1orig=R1;
clear R1;
R1.Chr=[];
R1.Start=[];
R1.End=[];
l=0;
for j=1:length(chr)
    chr_idx=strcmp(L.Chr{chr(j)},R1orig.Chr);
    R1.Chr=[R1.Chr; R1orig.Chr(chr_idx)];
    R1.Start=[R1.Start; R1orig.Start(chr_idx)+l];
    R1.End=[R1.End; R1orig.End(chr_idx)+l];
    l=l+L.Len(chr(j));
end
clear R1orig;

R2orig=R2;
clear R2;
R2.Chr=[];
R2.Start=[];
R2.End=[];
l=0;
for j=1:length(chr)
    chr_idx=strcmp(L.Chr{chr(j)},R2orig.Chr);
    R2.Chr=[R2.Chr; R2orig.Chr(chr_idx)];
    R2.Start=[R2.Start; R2orig.Start(chr_idx)+l];
    R2.End=[R2.End; R2orig.End(chr_idx)+l];
    l=l+L.Len(chr(j));
end
clear R2orig;

xmax=l;
Map.Chr=L.Chr(chr);
Map.Start=1;
Map.End=xmax;

[~,R1]=preprocessR(Map,R1);
[~,R2]=preprocessR(Map,R2);

card1=sum(R1{1}(:,2)-R1{1}(:,1)+1);
card2=sum(R2{1}(:,2)-R2{1}(:,1)+1);

nR1=length(R1{1}(:,1));
nR2=length(R2{1}(:,1));
r1(1:2:2*nR1-1)=R1{1}(:,1);
r1(2:2:2*nR1)=R1{1}(:,2);
r2(1:2:2*nR2-1)=R2{1}(:,1);
r2(2:2:2*nR2)=R2{1}(:,2);
r12=range_intersection(r1,r2);
card12=sum(r12(2:2:end)-r12(1:2:end-1)+1);


