clear;
close all;
clc;



%% Select chromosome and reference sets
addpath('C:\Users\Darshan\Documents\MATLAB\Transposable Elements\TE Insertion Bias\Data')
load Chr_lengths_human
% load mappability_human
load reference_human_cancer_new_comparison4
Rfull=R; 
clear R;

% Select reference sets
R(1)=Rfull(4);
R(2)=Rfull(5);
R(3)=Rfull(7);
R(4)=Rfull(8);
R(5)=Rfull(9);

chr=19; % Try 18 too
% chr=1:24;

N=length(R);
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

xmax=l;
% x=1:xmax;

% MapOrig=Map; 
% clear Map;
% l=0;
% Map.Chr=[];
% Map.Start=[];
% Map.End=[];
% for i=1:length(chr)
%     chr_idx=strcmp(L.Chr{chr(i)},MapOrig.Chr);
%     Map.Chr=[Map.Chr; MapOrig.Chr(chr_idx)];
%     Map.Start=[Map.Start; MapOrig.Start(chr_idx)+l];
%     Map.End=[Map.End; MapOrig.End(chr_idx)+l];
%     l=l+L.Len(chr(i));
% end
% clear MapOrig;

Map.Chr=L.Chr(chr);
Map.Start=1;
Map.End=xmax;

[Map,R]=preprocessR(Map,R);
M=formMstruct(R,Map);

%% Create pi (set indexing scheme), set theta, compute Q and g

pow2=2.^(0:N-1);
combos=dec2bin(1:2^N-1);
pi=zeros(2^N-1,N);
npi=zeros(2^N-1,1);
for i=1:2^N-1       
    for j=1:N
        pi(i,j)=str2double(combos(i,j));
    end
    npi(i)=pow2*pi(i,:)';
end
pi=pi(npi,:);

theta=zeros(2^N-1,1);
pi_true=[0 0 1 0 0;...
         0 0 0 1 0;...
         0 0 1 1 0];
theta(pow2*pi_true(1,:)')=1;
theta(pow2*pi_true(2,:)')=0.5;
theta(pow2*pi_true(3,:)')=-1;
param_ind_true=zeros(2^N-1,1);
param_ind_true(theta~=0)=1;
D=length(pi_true(:,1));

Q=createQ(N);
g=Q*theta;


%% Simulate query intervals and estimate paramters theta_hat

numIter=1000; % Number of replicates
n=2540; % Number of query intervals simulated per replicate
% n=10000;

theta_hat=zeros(D,numIter);
piout{numIter}=[];
% param_ind=ones(2^N-1,1);
% sum_param_ind=zeros(2^N-1,1);
h=waitbar(0,'Simulating Data & Selecting Model');
tic

%%%% This goes inside the for-loop if the distribution of y is not ones%%%%
%
% Set query interval lengths
qy=ones(n,1); % Change this to sample from some distribution if desired
[m,y]=hist(qy,unique(qy));
y(m==0)=[];
m(m==0)=[];
d=length(y);
%
% Compute S
S=computeS_f0uniform(M,y);
%
% Compute c
c=1./(S*exp(Q*theta));
%
% Null distribution (assuming all x is mappable)
f0=zeros(d,1);
for i=1:d
    f0(i)=1/(Map(2)-y(i)+1);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
for ii=1:numIter
    
    disp(num2str(ii));
    
    % Sample qx by acceptance/rejection
    qx=zeros(n,1);
    rq=zeros(n,2^N-1);
    for i=1:n

        yidx=find(y==qy(i));
        f0y=f0(yidx);
        cy=c(yidx);
        f=f0y*cy*exp(g);
        p=f/max(f);

        flag=true;
        while flag
            u=rand;
            xtry=round(xmax*rand+0.5);
            b=computeb(xtry,qy(i),R);
            if u<p(pow2*b'+1)
                qx(i)=xtry;
                rq(i,:)=getRowQ(b,Q);
                flag=false;
            end
        end

    end

    % Compute T
    T=sum(rq);
    
    % Compute MLE of true model 
    [theta_hat_full,~]=computeMLE(param_ind_true,m,T,S,Q);
    theta_hat(:,ii)=theta_hat_full(param_ind_true==1);
    
    % Do model selection
    [~,param_ind]=stepwise('SBC',ones(2^N-1,1),m,T,S,Q);
    piout{ii}=pi(param_ind==1,:);
%     sum_param_ind=sum_param_ind+param_ind;
    % pvals=computePvals(param_ind,m,T,S,Q);
    % CI=computeCI(theta_hat,param_ind,m,S,Q);
    
    waitbar(ii/numIter);
    
end
toc
close(h);

% Histograms of all model terms
nbins=round(sqrt(numIter));
for i=1:D
    figure(i); clf; hold on;
    hist(theta_hat(i,:),nbins);
    disp(['Mean of param ' num2str(i) ' is ' num2str(mean(theta_hat(i,:)))]);
    disp(['Std of param ' num2str(i) ' is ' num2str(std(theta_hat(i,:)))]);
    ylims=get(gca,'Ylim');
    plot([theta(pow2*pi_true(i,:)') theta(pow2*pi_true(i,:)')],ylims,'r','LineWidth',2);
end

% % Bar graph of percentage of times a model term is included
% figure(2^N); bar(sum_param_ind/numIter*100);
% xlim([0 2^N])

% Histogram or abundance plot of piout (this is going to be tough)
modelstr{numIter}=[];
for ii=1:numIter
    [M,~]=size(piout{ii});
    for m=1:M
        idx=find(piout{ii}(m,:)==1);
        s='{';
        for j=1:length(idx)
            s=[s num2str(idx(j)) ','];
        end
        s(end)='}';
        modelstr{ii}=[modelstr{ii} s ','];
    end
    modelstr{ii}(end)=[];
end
[C,~,ic]=unique(modelstr);


[idx,label]=grp2idx(sort(modelstr));
h=hist(idx,unique(idx));
% set(gca,'xTickLabel',label); % Need to rotate x-labels
modelperc=h/numIter*100;
[modelperc_sort,sortidx]=sort(modelperc,'descend');
csum=cumsum(modelperc_sort);
idx=find(csum>95);
Ctop95=C(sortidx(1:idx(1)));
modelperc_sort_top95=modelperc_sort(1:idx(1));
for i=1:idx(1)
    disp([Ctop95{i} ': ' num2str(modelperc_sort_top95(i))]);
end

% Save results
% save('simulation2_1000reps.mat');
