%% SCRIPT_TestArcLengthParamPP
clear all
close all
clc

%% Setup fit
n = 50;
% x = linspace(0,10,n);
% Y = rand(2,n);
x = linspace(-pi,pi,n);
Y(1,:) = 5*cos(x);
Y(2,:) = 10*sin(x);

pp = spline(x,Y);

%% Plot fit
N = 10000;
xx = linspace(min(x),max(x),N);
YY = ppval(pp,xx);

% -> Check tangent
ppT = diffpp(pp);
T = ppval(ppT,xx);
normT = sqrt( sum(T.^2,1) );

figure;
plot(xx,normT,'.b');

% -> Create 2D figure
fig2D(1) = figure;
axs2D(1) = axes('Parent',fig2D(1));
hold(axs2D(1),'on');
xlabel(axs2D(1),'x');
ylabel(axs2D(1),'y');

% -> Plot 2D fit
plt2D(1) = plot(axs2D(1),YY(1,:),YY(2,:),'b');

% -> Plot fit evolving as a function of x
fig(2) = figure;
m = size(Y,1);
ylbls = {'x','y'};
for i = 1:m
    axs(1,i) = subplot(m,1,i);
    hold(axs(1,i),'on');
    ylabel(axs(1,i),ylbls{i})
    
    % -> Plot "time" evolving variable
    plt(i) = plot(axs(1,i),xx,YY(i,:),'b');
end

%% Parametrize for arc length
[pps,slim] = arcLengthParamPP(pp,50);
ss = linspace(slim(1),slim(2),N);
YYs = ppval(pps,ss);

% -> Plot 2D fit
plt2D(1) = plot(axs2D(1),YYs(1,:),YYs(2,:),'r');

% -> Plot fit evolving as a function of s
for i = 1:m
    % -> Plot "time" evolving variable
    plt(i) = plot(axs(1,i),ss,YYs(i,:),'r');
end

% -> Check tangent
ppTs = diffpp(pps);
Ts = ppval(ppTs,ss);
normTs = sqrt( sum(Ts.^2,1) );

figure;
plot(ss,normTs,'.r');