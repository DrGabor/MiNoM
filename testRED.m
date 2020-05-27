clc; close all; clear all; 
Err = importdata('residuals/1.txt');
P = [1.0 2.0]; 
maxIter_em = 10; 
bRand = 1; 
[Model, Z] = MoEPFittingFun(Err, P, maxIter_em, bRand);
Edges = linspace(0, max(Err), 50); 
NoiseAnalysisFun(Model, Err, Edges); 