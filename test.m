% Script to test the MatrixSOPOT class.


clc
clear
close all


a = MatrixSOPOT([.345+1j*.3 .456+1j*.6 .567+1j*.4; .123 .234 .345], 4, 16, 3);
b = MatrixSOPOT(.55, 4, 16, 3);
% EoF
