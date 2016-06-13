clear;
close all;
clc;

% Convert excel sheet (bed file) to mat file

%% Create query or reference interval set

query_or_reference='reference';
output_path='C:\Users\Darshan\Desktop\GINOM\mat_files';

if strcmp(query_or_reference,'query')
    
    data_path='C:\Users\Darshan\Desktop\GINOM\raw_data';
    input_filename='query_som_trans_lung.xlsx';
    output_filename='query_ex.mat';
    Q = CreateQueryMatFromExcel(input_filename,output_filename,data_path,output_path);
    
elseif strcmp(query_or_reference,'reference')
    
    data_path='C:\Users\Darshan\Desktop\GINOM\raw_data';
    input_filename='Broad_A549_H3k04me3_Dex.xlsx,reference_2_hg19_Refseq_promoters_2000bp.xlsx,reference_3_hg19_Refseq_gene_bodies.xlsx,Broad_A549_H3k36me3_Dex.xlsx,Broad_A549_H3k79me2_Dex.xlsx';
    output_filename='reference_ex.mat';
    R = CreateReferenceMatFromExcel(input_filename,output_filename,data_path,output_path);
    
end

%% Create mappability mat file and Chr length file


