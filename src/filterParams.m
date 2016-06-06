function param_ind = filterParams(pi,maxorder,inputcombos)

[~,N]=size(pi);

if ~strcmp(maxorder,'NA')
    mult=sum(pi,2);
    param_ind_1=mult<=str2double(maxorder);
else
    param_ind_1=true(2^N-1,1);
end

if ~strcmp(inputcombos,'NA')
    [a,~]=size(inputcombos);
    pow2=2.^(0:N-1);
    combos_num=zeros(a,N);
    npi_input=zeros(a,1);
    for i=1:a        % (there's probably an efficient way to do this)
        for j=1:N
            combos_num(i,j)=str2double(inputcombos(i,j));
        end
        npi_input(i)=pow2*combos_num(i,:)';
    end
    param_ind_2=zeros(2^N-1,1);
    param_ind_2(npi_input)=1;
    param_ind_2=logical(param_ind_2);
else
    param_ind_2=true(2^N-1,1);
end

param_ind=double(param_ind_1 & param_ind_2);

