# GINOM (Genomic INterval Overlap Model)
GINOM implements a statistical framework for assessing interval overlap of multiple genomic features

Example workflow:

Step 1) Convert .bed files to .mat files. 

Download a .bed file for the desired query interval set as well as a .bed file for each reference interval set. Save each .bed file as an excel file (.xlsx). Use the script run_CreateMatFromExcel.m to convert these spreadsheets to .mat files, which will be used later as input to a data preprocessing function. Instructions for run_CreateMatFromExcel.m are the following: 

a) Set query_or_reference equal to 'query', and then set output_path equal to the path where you wish to store the query interval .mat file to be created. Within the first if-statement, set data_path equal to the path where your query interval excel spreadsheet is stored, and set input_filename to the query interval file name in the format 'filename.xlsx'. Set output_filename to the file name of the .mat file to be created in the format 'filename.mat'. Run the script. The function CreateQueryMatFromExcel.m will perform the file conversion and save the .mat file in the specified location. 

The .mat file will hold a structure array Q with the following fields: Chr, Start, and End. Start and End are both n x 1 arrays that contain the starting and ending nucleotide index of each of the n query intervals, respectively. Chr is an n x 1 cell array that indicates on which chromosome each of the n query intervals is located.    

b) Set query_or_reference equal to 'reference', and then set output_path equal to the path where you wish to store the reference interval .mat file to be created. Within the elseif-statement, set data_path equal to the path where your reference interval excel spreadsheets are stored. Set input_filename as a string that consists of all reference interval file names, each separated by commas (no spaces) as in for example 'filename1.xlsx,filename2.xlsx,filename3.xlsx'. The order of the filenames here sets the indexing order of the reference sets from here on. Set output_filename to the file name of the .mat file to be created in the format 'filename.mat'. Run the script. The function CreateQueryMatFromExcel.m will perform the file conversion and save the .mat file in the specified location.  

The .mat file will hold a structure R that consists of N structure arrays, each corresponding to a reference set in the order specified above. Structure R(i), the ith reference set, has the same fields as described in Q. 

Notice that there is an empty space in the run_CreateMatFromExcel.m script allocated for creating the mappability file and the chromosome length file. For now, we have provided these files for the human genome (mm9), and in a later update we will fill out the rest of this script to automate the creation of these .mat files. For now, include in your data folder the following two files: mappability_human.mat and Chr_lengths_human.mat. The mappability file holds a structure called Map, which consists of the same fields as the reference and query intervals as explained above. The intervals in this case represent the intervals of the genome that have non-zero mappability. The file Chr_lengths_human.mat holds a structure L with two fields, Chr and Len. Len is a 25 x 1 array consisting of the length, or the number of nucleotides, of chromosomes 1-22, X, Y, and M, in that order. Chr is a 25 x 1 cell array that labels each index of Len with the appropriate chromosome name. 

We now have all the raw data files necessary to run preprocessing.



Step 2) Preprocess the raw data 

The data preprocessing is handled by the function preprocessData.m, which takes as input a .txt file that consists of specially formatted instructions. This configuration file consists of 8 lines of text, and the contents of each line is explained in the file preprocess_legend.txt. We have also provided an example configuration file, preprocess_NC1.txt, to help explain the formatting. Here, we provide more detailed instruction on how to create this configuration file for the preprocessing step.  

Line 1: The file name of the reference interval set, as was output from running the script run_CreateMatFromExcel.m. Ex: reference.mat

Line 2: The file name of the query interval set, as was output from running the script run_CreateMatFromExcel.m. Ex: query.mat

Line 3: The file name of the mappable intervals. We provide this file. Ex: mappability_human.mat

Line 4: The file name of the chromosome lengths. We provide this file. Ex: Chr_lengths_human.mat

Line 5: The path containing the above .mat data files. Ex: C:\Users\John\Documents\MATLAB\Data

Line 6: The path to save the output data file consisting of the preprocessed data. Ex: C:\Users\John\Documents\MATLAB\Data

Line 7: The name of the output data file consisting of the preprocessed data. Ex: preprocessedData.mat

Line 8: The indices of chromosomes that one wishes to include in the analysis. The indexing follows the labeling convention provided in L.Chr from the chromosome lengths data file (Chr_lengths_human.mat). The format here is to list all indices as numbers separated by commas (no spaces). Ex: 1,2,3,4,5,6 to analyze chromosomes 1-6.    

After creating the preprocessing configuration text file as shown above, save it in your working Matlab directory. To run the preprocessing only, open the script run_model.m, uncomment the preprocessData() function call, and comment the fitModel() function call. Type your preprocessing configuration text file name in single quotes as input to the preprocessData() function. After the preprocessing function finishes running, look for the .mat file to have appeared in your data folder as specified by lines 6 and 7 of the configuration text file. This .mat file consists of the following variables: 

pi: Set indexing scheme, notated by \pi in Section 2.1 of the main text

b: Binary indexing scheme, notated by \beta in Section 2.4 of the main text

m: Multiplicity of query interval lengths, as described in supplemental material

Q: Described in the supplemental material

S: Described in the supplemental material

T: Vector of overlap counts, indexed by pi and described above Eqn. 4 of the main text

The preprocessing step is complete, and now we can move on to fitting the model. 



Step 3) Fit the model


The model fitting is handled by the function fitModel.m, which takes as input a .txt file that consists of specially formatted instructions, as in the preprocessing case. This configuration file consists of 7 lines of text, and the contents of each line is explained in the file fitmodel_legend.txt. We have also provided an example configuration file, fitmodel_NC1.txt, to help explain the formatting. Here, we provide more detailed instruction on how to create this configuration file for the model fitting step.  

Line 1: The file name of the preprocessed data, as output from Step 2. This line should be the same as line 7 of the preprocessing configuration text file. Ex: preprocessedData.mat

Line 2: The path containing the .mat data files. Ex: C:\Users\John\Documents\MATLAB\Data 

Line 3: The path to save the fitModel.m output file. Ex: C:\Users\John\Documents\MATLAB\Results

Line 4: The output file name. Ex: modelfit.mat

Line 5: Select the information criterion you wish to use for stepwise model selection (see Section 2.3 of the main text). There are three options for this line: SBC, AIC, or NA. The SBC penalty is more selective than the AIC penalty, and the fitted model will generally contain less model terms. Set to NA if you do not wish to do any model selection. EX: SBC

Line 6: The maximum order of interaction term you wish to include in the fitted model. This entry will be a single number, or it will be NA if you do not wish to impose any such restriction. Ex: 3. In this case, \pi will be restricted to cardinality 3 or less in all terms \theta_{\pi} to be considered in the model.

Line 7: Explicitly select the model terms you wish to consider for model fitting. Set to NA to not impose any such restriction; otherwise, list the allowed model terms encoded as binary strings, one on each line. Ex: 
10000
01000
00100
00010
00001
01100
00011
This sequence indicates that the model terms \theta_{\pi} with \pi equal to {1}, {2}, {3}, {4}, {5}, {2,3}, and {4,5} are allowed for consideration in the model fitting procedure.  

Here, we provide some examples on how to use lines 5, 6, and 7 to configure the options in fitModel.m. 

a) Use stepwise model selection with SBC penalty to fit a model consisting of only 3rd order terms or less: 

Line 5: SBC
Line 6: 3
Line 7: NA

b) Use stepwise model selection with SBC penalty to fit a model with no other restrictions

Line 5: SBC
Line 6: NA
Line 7: NA

c) Fit the full model:

Line 5: NA
Line 6: NA
Line 7: NA

d) Fit the model that consists of all first order terms: 

Line 5: NA
Line 6: 1
Line 7: NA

or (assuming 5 reference sets)

Line 5: NA
Line 6: NA
Line 7: 
10000
01000
00100
00010
00001

e) Fit the model that consists of only terms {1} and {2,3} (assuming 5 reference sets):

Line 5: NA
Line 6: NA
Line 7: 
10000
01100

f) Use stepwise model selection with SBC penalty to select the best model from the terms {1}, {2}, {3}, {4}, {1,2}, {2,3}, and {3,4}:

Line 5: SBC
Line 6: NA
Line 7: 
10000
01000
00100
00010
11000
01100
00110

After creating the fitmodel configuration text file as shown above, save it in your working Matlab directory. To fit the model, open the run_model.m script, comment the line with the preprocessData function call, uncomment the line with the fitModel function call, enter the name of your fitmodel configuration text file in single quotes as input to fitModel(), and run the script. Of course, you can run both the preprocessing and the model fitting by uncommenting both lines of the run_model.m script. The output of fitModel() consists of 6 variables. 

pi: The set indexing scheme of the model parameters included in the fitted model, notated by \pi in Section 2.1 of the main text. Each row indicates a model parameter \pi. 

theta_hat: The value of the fitted model parameters, indexed according to the rows of pi. 

pvals: The p-values associated with each model parameter, indexed by pi.

CI: The 95% confidence interval of each value of theta_hat. 

b: Binary indexing scheme, notated by \beta in Section 2.4 of the main text. 

g: The enrichment function, described in Section 2.4, indexed by the rows of b. 

This concludes the tutorial on GINOM. Now, we present a brief set of instructions on how to run the example provided. 



Running the Example:

Open the run_CreateMatFromExcel.m script. Change the output_path variable to reflect the folder on your machine where you wish for the .mat files to be stored. This folder should also contain mappability_human.mat and Chr_lengths_human.mat. Change the data_path variable to reflect the folder on your machine where the excel files are stored. The input and output file names are already set for this example. Set query_or_reference equal to 'query', and run the script. Set query_or_reference equal to 'reference', and run the script again. 

Open preprocess_ex.txt. Change line 5 to reflect the folder on your machine where the .mat data files are stored (this will be the same as output_path in the run_CreateMatFromExcel.m script). Change line 6 to reflect the folder on your machine where you wish the preprocessed data to be stored (it can be the same as line 5). Save preprocess_ex.txt to your working directory in Matlab.

Open fitmodel_ex.txt. Change line 2 to reflect the folder on your machine where the preprocessed data is stored. Change line 3 to reflect the folder on your machine where you wish for the fitted model results to be stored. Save fitmodel_ex.txt to your working directory in Matlab.

Run the script run_model.m. Total computation time will be around 25 seconds on a 2.5gHz laptop. The output variables pi and theta_hat will equal [0 0 0 1 0] and -0.8934, respectively. This means that only model term {4} is included in the model and that query intervals are exp(-0.8934)=0.4093 times more likely to overlap reference set 4 than what would be expected under the null hypothesis. 

If you wish to run on a different set of chromosomes, change line 8 of preprocess_ex.txt accordingly, save, and re-run the run_model.m script. If, for example, you wish to run on all chromosomes, change line 8 to 
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24. The result will be 
pi=[0 1 0 0 0; 
    0 0 0 1 0;
    0 0 0 0 1;
    0 0 0 1 1], 
and 
theta_hat=[-0.5683;
           -0.7135;
           -0.8756;
            0.6578]. 
This analysis will take approximately 20 minutes to run in this case. The terms included in the model are {2}, {4}, {5}, and {4,5}. See Section 2.4 of the main text for more information on model interpretation. 


