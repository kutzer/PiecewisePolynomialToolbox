%% SCRIPT_TestPPfunctions
clear all
close all
clc

%% Define candidate functions and samplingand plot candidate data
% -> Define loose sampling
n = 50;
x = linspace(-pi/3,4*(2*pi),n);
% -> Define tight sampling
nn = 20*n;
xx = linspace(min(x),max(x),20*n);

% -> Define functions to evaluate
f{1}  = @(x_in)  4*sin(x_in) + 10;
df{1} = @(x_in)  4*cos(x_in);

f{2}  = @(x_in)  3*cos(x_in) + 8;
df{2} = @(x_in) -3*sin(x_in);

f{3}  = @(x_in) 1*sin(x_in).^2 + 6;
df{3} = @(x_in) 2*sin(x_in).*cos(x_in);

%% Create figures and axes
ylbls = {'f(x)','\frac{df(x)}{dx}'};
for i = 1:numel(f)
    fig(i) = figure('Name',sprintf('Function %d: %s',i,char(f{i})));
    for j = 1:2
        axs(j,i) = subplot(2,1,j,'Parent',fig(i));
        hold(axs(j,i),'on');
        xlabel(axs(j,i),'x');
        ylabel(axs(j,i),ylbls{j});
    end
end

% -> Define output data for each function
for i = 1:numel(f)
    y(i,:)  = f{i}(x);
    yy(i,:) = f{i}(xx);
end

% -> Plot output data for each function
for i = 1:size(y,1)
    pltY(i) = plot(axs(1,i),x,y(i,:),'xr');
    pltYY(i) = plot(axs(1,i),xx,yy(i,:),'r');
end

%% Define and plot derivative
% -> Define output data for each function derivative
for i = 1:numel(df)
    dy(i,:) = df{i}(x);
    dyy(i,:) = df{i}(xx);
end

% -> Plot output data for each function derivative
for i = 1:size(dy,1)
    pltDY(i) = plot(axs(2,i),x,dy(i,:),'xr');
    pltDYY(i) = plot(axs(2,i),xx,dyy(i,:),'r');
end

%% Define and plot piecewise polynomial fit
pp = spline(x,y);
yyPP = ppval(pp,xx);

for i = 1:size(yyPP,1)
    pltPP(i) = plot(axs(1,i),xx,yyPP(i,:),'b');
end

%% Define and plot derivative of piecewise polynomial fit
dpp = diffpp(pp);
dyyPP = ppval(dpp,xx);

for i = 1:size(dyyPP,1)
    pltDPP(i) = plot(axs(2,i),xx,dyyPP(i,:),'b');
end

%% Integrate piecewise polynomial derivative and plot
idpp = intpp(dpp);
idyyPP = ppval(idpp,xx);

for i = 1:size(idyyPP,1)
    pltIDPP(i) = plot(axs(1,i),xx,idyyPP(i,:),'m');
end

%% Apply integration constatnt at first breakpoint
idpp = intpp(dpp,y(:,1));
idyyPP = ppval(idpp,xx);

for i = 1:size(idyyPP,1)
    pltIDPP(i) = plot(axs(1,i),xx,idyyPP(i,:),'g:','LineWidth',1.7);
end

%% Apply integration constant at random breakpoint
x_rand = (max(x) - min(x)).*rand(1) + min(x);
y_rand = ppval(pp,x_rand);
idpp = intpp(dpp,y_rand,x_rand);
idyyPP = ppval(idpp,xx);

for i = 1:size(idyyPP,1)
    pltIDPP(i) = plot(axs(1,i),xx,idyyPP(i,:),'c--','LineWidth',1);
end