function [theta_hat,pvals,CI,g,pi,b] = fitModel(filename) %#ok<STOUT>

disp('Reading input text file');
tic
timeidx=0;

fid=fopen(filename);
C=textscan(fid,'%s','Delimiter','|');

data_name=C{1}{1};
data_path=C{1}{2};
output_path=C{1}{3};
output_name=C{1}{4};
selection_criterion=C{1}{5};
maxorder=C{1}{6};
L=length(C{1});
for i=7:L
    inputcombos(i-6,:)=C{1}{i}; %#ok<AGROW>
end

timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

disp('Loading data');
tic

% Add data path to working directory
addpath(data_path);

% Load reference interval sets, R 
load(data_name); 

% Filter params 
param_ind_init=filterParams(pi,maxorder,inputcombos); %#ok<NODEF>

timeidx=timeidx+1;
toc
time(timeidx)=toc;
disp(' ');

disp('Fitting model')
tic
% Fit model
if strcmp(selection_criterion,'NA')
    % No model selection. Fit the full model.
    [theta_hat,~]=computeMLE(param_ind_init,m,T,S,Q); 
    param_ind=param_ind_init;
else
    % Stepwise model selection.
    [theta_hat,param_ind]=stepwise(selection_criterion,param_ind_init,m,T,S,Q); 
    
    % Test for hidden zero-overlaps problem, and rerun if necessary
    [rerun,param_ind_rerun]=testForHidden0ovl(pi(param_ind==1,:),T,param_ind_init);
 
    if rerun
        timeidx=timeidx+1;
        toc
        time(timeidx)=toc;
        disp(' ');
        
        disp('Hidden zero overlap detected and removed. Rerunning model selection...');
        tic
        [theta_hat,param_ind]=stepwise(selection_criterion,param_ind_rerun,m,T,S,Q);
    end
    
%     % Try all model combinations
%     [theta_hat,param_ind]=modelSelectFromAll(selection_criterion,param_ind,m,T,S,Q);
end

%keyboard
if ~isempty(theta_hat)
    
    g=Q*theta_hat;
    timeidx=timeidx+1;
    toc
    time(timeidx)=toc;
    disp(' ');

    % Compute p-values of model parameters
    disp('Computing p-values')
    tic
    pvals=computePvals(param_ind,m,T,S,Q);
    timeidx=timeidx+1;
    toc
    time(timeidx)=toc;
    disp(' ');

    % Compute confidence intervals
    disp('Computing confidence intervals')
    tic
    CI=computeCI(theta_hat,param_ind,m,S,Q);
    timeidx=timeidx+1;
    toc
    time(timeidx)=toc;
    disp(' ');
else
    g=[];
    pvals=[];
    CI=[];
    timeidx=timeidx+1;
    toc
    time(timeidx)=toc;
end

% Save output
disp('Saving results')
tic
theta_hat=theta_hat(param_ind==1);
pi=pi(param_ind==1,:);
save([output_path output_name],'theta_hat','pvals','CI','g','pi','b');
toc
time(timeidx)=toc;
disp(' ');
disp(['Total elapsed time is ' num2str(sum(time)) ' seconds.']) 
