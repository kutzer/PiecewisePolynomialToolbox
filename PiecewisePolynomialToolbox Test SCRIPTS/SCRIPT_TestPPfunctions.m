%% SCRIPT_TestPPfunctions
clear all
close all
clc

%% Create figure and axes
fig = figure;
axs(1) = subplot(2,1,1,'Parent',fig);
axs(2) = subplot(2,1,2,'Parent',fig);
hold(axs(1),'on');
hold(axs(2),'on');

%% Define and plot candidate data
n = 50;
x = linspace(-pi/3,4*(2*pi),n);
xx = linspace(min(x),max(x),20*n);

y = 4*sin(x) + 10;
yy = 4*sin(xx) + 10;

pltY = plot(axs(1),x,y,'xr');
pltYY = plot(axs(1),xx,yy,'r');

%% Define and plot derivative
dy = 4*cos(x);
dyy = 4*cos(xx);

pltDY = plot(axs(2),x,dy,'xr');
pltDYY = plot(axs(2),xx,dyy,'r');

%% Define and plot piecewise polynomial fit
pp = spline(x,y);
yyPP = ppval(pp,xx);

pltPP = plot(axs(1),xx,yyPP,'b');

%% Define and plot derivative of piecewise polynomial fit
dpp = diffpp(pp);
dyyPP = ppval(dpp,xx);

pltDPP = plot(axs(2),xx,dyyPP,'b');

%% Integrate piecewise polynomial derivative and plot
idpp = intpp(dpp);
idyyPP = ppval(idpp,xx);

pltIDPP = plot(axs(1),xx,idyyPP,'m');

%% Apply integration constatnt at first breakpoint
idpp = intpp(dpp,y(1));
idyyPP = ppval(idpp,xx);

pltIDPP = plot(axs(1),xx,idyyPP,'g:','LineWidth',1.7);

%% Apply integration constant at random breakpoint
x_rand = (max(x) - min(x)).*rand(1) + min(x);
y_rand = ppval(pp,x_rand);
idpp = intpp(dpp,y_rand,x_rand);
idyyPP = ppval(idpp,xx);

pltIDPP = plot(axs(1),xx,idyyPP,'c--','LineWidth',1);