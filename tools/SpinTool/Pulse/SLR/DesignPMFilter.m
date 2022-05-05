clc
clear all
close all


B=firceqrip(128, 0.4, [0.001 0.001]); [h,w] = freqz(B,1,501); figure; plot(abs(h));
