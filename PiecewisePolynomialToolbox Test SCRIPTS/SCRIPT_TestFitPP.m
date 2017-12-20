%% SCRIPT_TestFitPP
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

%% FitPP (no constraints on derivative)
xxx = linspace(min(x),max(x),50000);
fpp = fitpp(x,y);
yyFPP = ppval(fpp,xxx);

pltFPP = plot(axs(1),xxx,yyFPP,'k');

%% Define and plot derivative of piecewise polynomial fit
dfpp = diffpp(fpp);
dyyFPP = ppval(dfpp,xxx);

pltDFPP = plot(axs(2),xxx,dyyFPP,'b');

%% Integrate piecewise polynomial derivative and plot
idfpp = intpp(dfpp);
idyyFPP = ppval(idfpp,xxx);

pltIDFPP = plot(axs(1),xxx,idyyFPP,'m');

%% Apply integration constatnt at first breakpoint
idfpp = intpp(dfpp,y(1));
idyyFPP = ppval(idfpp,xxx);

pltIDFPP = plot(axs(1),xxx,idyyFPP,'g--','LineWidth',1.5);

%% Apply integration constant at random breakpoint
x_rand = (max(x) - min(x)).*rand(1) + min(x);
y_rand = ppval(fpp,x_rand);
idfpp = intpp(dfpp,y_rand,x_rand);
idyyFPP = ppval(idfpp,xxx);

pltIDFPP = plot(axs(1),xxx,idyyFPP,'c--','LineWidth',1);

%% FitPP (with constraints on derivative)
fpp = fitpp(x,y,x,dy);
yyFPP = ppval(fpp,xx);

pltFPP = plot(axs(1),xx,yyFPP,'b');

%% Define and plot derivative of piecewise polynomial fit
dfpp = diffpp(fpp);
dyyFPP = ppval(dfpp,xx);

pltDFPP = plot(axs(2),xx,dyyFPP,'b');

%% Integrate piecewise polynomial derivative and plot
idfpp = intpp(dfpp);
idyyFPP = ppval(idfpp,xx);

pltIDFPP = plot(axs(1),xx,idyyFPP,'m');

%% Apply integration constatnt at first breakpoint
idfpp = intpp(dfpp,y(1));
idyyFPP = ppval(idfpp,xx);

pltIDFPP = plot(axs(1),xx,idyyFPP,'g--','LineWidth',1.5);

%% Apply integration constant at random breakpoint
x_rand = (max(x) - min(x)).*rand(1) + min(x);
y_rand = ppval(fpp,x_rand);
idfpp = intpp(dfpp,y_rand,x_rand);

pltIDFPP = plot(axs(1),xx,idyyFPP,'c--','LineWidth',1);