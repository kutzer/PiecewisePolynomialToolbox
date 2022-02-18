%% SCRIPT_Test_ppArray2pp
clear all
close all
clc

%% Create piecewise polynomials
% Define points & independent variable for 1D piecewise fit
x1 = rand(1,5);                   % 1D points
s1 = linspace(0,1,size(x1,2));    % Independent variable
% Fit piecewise polynomial (cubic spline)
ppArray(1) = spline(s1,x1);

% Define points & independent variable for 3D piecewise fit
x2 = rand(3,7);                   % 3D points
s2 = linspace(-1,2,size(x2,2));   % Independent variable
% Fit piecewise polynomial (C^0 continuous)
ppArray(2) = fitpp(s2,x2);

% Define points & independent variable for 3D piecewise fit
x3 = rand(3,10);                   % 3D points
dx3 = rand(3,10);    
ddx3 = rand(3,10);
s3 = linspace(-0.5,2,size(x3,2));   % Independent variable
% Fit piecewise polynomial (C^0 continuous)
ppArray(3) = fitpp(s3,x3,s3,dx3,s3,ddx3);

%% Combine piecewise polynomial array to a single piecewise polynomial
pp = ppArray2pp(ppArray);