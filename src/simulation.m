clear;
close all;
clc;

xmax=100;
x=1:xmax;
Map=[1 xmax];


R{1}=[21 60];
% R{2}=R{1};
R{2}=[41 80];
% R{2}=[41 80];
% R{1}=[11,30;41,60;71,90];
% R{2}=[6,25;33,38;66,85;93,98];
% R{3}=[16,35;76,95];
N=length(R);


% theta=[0.5;0.7;-0.2];
theta=[-20;1;21];
% theta=[1;0;0;1;0;0;1];

n=2540;

qy=ones(n,1); % Change this to sample from some distribution if desired

[m,y]=hist(qy,unique(qy));
y(m==0)=[];
m(m==0)=[];
d=length(y);



f0=zeros(d,length(x));
for i=1:d
    for j=1:length(Map(:,1))
        if Map(j,2)-y(i)+1>=Map(j,1)
            f0(i,Map(j,1):Map(j,2)-y(i)+1)=1;
        end
    end
    f0(i,:)=f0(i,:)/sum(f0(i,:));
end

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

b{d}=[];
for i=1:d
    b{i}=zeros(N,length(x));
    for j=1:length(x)
        b{i}(:,j)=computeb(x(j),y(i),R);
    end
end
r{d}=[];
for i=1:d
    r{i}=ones(2^N-1,length(x));
    for j=1:2^N-1
        ind=find(pi(j,:)==1);
        for k=1:length(ind)
            r{i}(j,:)=b{i}(ind(k),:).*r{i}(j,:);
        end
    end
end

numIter=1;
theta_hat=zeros(2^N-1,numIter);
for ii=1:numIter
    disp(num2str(ii));
    qx=zeros(n,1);
    for i=1:n

        yidx=find(y==qy(i));
        f0y=f0(yidx,:);
        ry=r{yidx};
        theta_ry=zeros(1,length(x));
        for k=1:2^N-1
            theta_ry=theta_ry+theta(k)*ry(k,:);
        end
        f=f0y.*exp(theta_ry);
        f=f/sum(f);
        p=f/max(f);

        flag=true;
        while flag
            u=rand;
            xtry=round(xmax*rand+0.5);
            if u<p(xtry)
                qx(i)=xtry;
                flag=false;
            end
        end

    end

    figure(1); clf; hold on; 
    hist(qx,x);
%     plot(x,n*f0,'g.','MarkerSize',15);
    plot(x,n*f,'r.','MarkerSize',15);
    plot([21,60],[-8,-8],'c','LineWidth',12);
    plot([41,80],[-8,-8],'m','LineWidth',6);
    set(gca,'FontSize',20,'FontName','Times');
    text(29,-15,'R_1','FontSize',20,'FontName','Times');
    text(69,-15,'R_2','FontSize',20,'FontName','Times');
    xlabel('x','FontSize',24,'FontName','Times'); 
    xlim([0,xmax+1]);

    Q=createQ(N);
    rq=zeros(n,2^N-1);
    for i=1:n
        rq(i,:)=getRowQ(computeb(qx(i),qy(i),R),Q);
    end
    T=sum(rq);
    M=formMstruct(R,Map);
    S=computeS_f0uniform(M,y);
    param_ind=ones(2^N-1,1);
%     [theta_hat(:,ii),~]=computeMLE(param_ind,m,T,S,Q);
    [theta_hat(:,ii),param_ind]=stepwise('SBC',param_ind,m,T,S,Q);
end

% pvals=computePvals(param_ind,m,T,S,Q);
% CI=computeCI(theta_hat,param_ind,m,S,Q);

% nbins=round(sqrt(numIter));
% 
% figure(2); clf; hold on; 
% hist(theta_hat(1,:),nbins);
% plot([theta(1) theta(1)],[0 50],'r','LineWidth',2);
% set(gca,'FontSize',16,'FontName','Times');
% xlim([0.4,0.6]);
% ylim([0,40]);
% 
% figure(3); clf; hold on; 
% hist(theta_hat(2,:),nbins);
% plot([theta(2) theta(2)],[0 50],'r','LineWidth',2);
% set(gca,'FontSize',16,'FontName','Times');
% xlim([0.6,0.8]);
% ylim([0,40]);
% 
% figure(4); clf; hold on; 
% hist(theta_hat(3,:),nbins);
% plot([theta(3) theta(3)],[0 50],'r','LineWidth',2);
% set(gca,'FontSize',16,'FontName','Times');
% xlim([-0.325,-0.075]);
% ylim([0,40]);


