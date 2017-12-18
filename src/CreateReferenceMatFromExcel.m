function R = CreateReferenceMatFromExcel(input_filename,output_filename,data_path,output_path)

addpath(data_path)

% Parse string
idx=strfind(input_filename,',');
if ~isempty(idx)
    idx=[0 idx length(input_filename)+1];
    n=length(idx)-1;
    ref_name{n}=[];
    for i=1:n
        ref_name{i}=input_filename(idx(i)+1:idx(i+1)-1);
    end
else
    ref_name{1}=input_filename;
end

% Create reference structure array
R=struct([]);
for i=1:n
    C=importdata(ref_name{i});
    R(i).Chr=C.textdata(:,1);
    R(i).Start=C.data(:,1);
    R(i).End=C.data(:,2);
end

save([output_path filesep output_filename],'R');