function Q = CreateQueryMatFromExcel(input_filename,output_filename,data_path,output_path)

addpath(data_path);
C=importdata(input_filename);
Q.Chr=C.textdata(:,1);
Q.Start=C.data(:,1);
Q.End=C.data(:,2);
save([output_path filesep output_filename],'Q');