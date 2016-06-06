clear;
close all;
clc;

% Convert excel sheet (bed file) to mat file

%% Create query or reference interval set

query_or_reference='reference';
output_path='C:\Users\Darshan\Documents\MATLAB\Transposable Elements\TE Insertion Bias\Data';

if strcmp(query_or_reference,'query')
    
    data_path='C:\Users\Darshan\Desktop\Work Research\Transposable Elements\TE Insertion Bias\Data\Human\Cancer\new_comparison1';
    input_filename='query_som_trans_lung.xlsx';
    output_filename='query_som_trans_lung.mat';
    Q = CreateQueryMatFromExcel(input_filename,output_filename,data_path,output_path);
    
elseif strcmp(query_or_reference,'reference')
    
    data_path='C:\Users\Darshan\Desktop\Work Research\Transposable Elements\TE Insertion Bias\Data\Human\Cancer\comparison_2';
    input_filename='wgEncodeAwgDnaseUwCaco2UniPk.xlsx,reference_2_hg19_Refseq_promoters_2000bp.xlsx,reference_3_hg19_Refseq_gene_bodies.xlsx,reference_4_hg19_pwmscan_L1HS.xlsx,reference_5_HCT116_H3K4me3.xlsx';
    output_filename='reference_human_cancer_comparison_2.mat';
    R = CreateReferenceMatFromExcel(input_filename,output_filename,data_path,output_path);
    
end

%% Create mappability mat file and Chr length file


