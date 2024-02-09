clc;clear;close all
load('F:\Github\GitHub\Tire-Data\Round 8\Cornering Run Data\B1965run1.mat')
figure(3)
plot(ET,IA)
figure(4)
plot(ET,SA)

load('F:\Github\GitHub\Tire-Data\Round 8\Cornering Run Data\B1965run2.mat')
figure(7)
plot(ET,IA)
figure(8)
plot(ET,SA)

load('F:\Github\GitHub\Tire-Data\Round 8\Cornering Run Data\B1965run4.mat')
figure(11)
plot(IA(8714:11502),FX(8714:11502))
figure(12)
plot(ET,SA)