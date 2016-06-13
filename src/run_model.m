clear;
close all;
clc;

preprocessData('preprocess_ex.txt');
[theta_hat,pvals,CI,g,pi,b]=fitModel('fitmodel_ex.txt');
