clear;
close all;
clc;

% preprocessData('preprocess_singleRef.txt');
[theta_hat,pvals,CI,g,pi,b]=fitModel('fitmodel_NC1.txt');
