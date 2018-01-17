%% SCRIPT_TestTangentPP
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

%% Define figure and axes
fig = figure;
axs = axes('Parent',fig);
hold(axs,'on');
daspect(axs,[1 1 1]);
view(axs,3);
xlabel(axs,'x');
ylabel(axs,'y');
zlabel(axs,'z');

%% Define spline
for i = 1:numel(f)
    Y(i,:) = f{i}(x);
    YY(i,:) = f{i}(xx);
end

pltYY = plot3(YY(1,:),YY(2,:),YY(3,:),'m');

pp = spline(x,Y);
ppY = ppval(pp,xx);
plotPP = plot3(ppY(1,:),ppY(2,:),ppY(3,:),'k');

%% Calculate tangent and unit tangent
[ppT, ppT_hat] = tangentpp(pp);

T_hat = ppval(ppT_hat,x);
for i = 1:numel(x)
    pltT_hat(i) = plot3(axs,[Y(1,i),Y(1,i)+T_hat(1,i)],[Y(2,i),Y(2,i)+T_hat(2,i)],[Y(3,i),Y(3,i)+T_hat(3,i)],'r');
end

%% Calculate magnitude of T_hat
TT_hat = ppval(ppT_hat,xx);
TT_hat_mag = sqrt(sum((TT_hat.^2),1));

figure;
plot(xx,TT_hat_mag);
xlabel('x');
ylabel('Magnitude, T hat');

%% Calculate normal and unit normal
[ppN, ppN_hat] = normalpp(pp,50,0);

N_hat = ppval(ppN_hat,x);
for i = 1:numel(x)
    pltT_hat(i) = plot3(axs,[Y(1,i),Y(1,i)+N_hat(1,i)],[Y(2,i),Y(2,i)+N_hat(2,i)],[Y(3,i),Y(3,i)+N_hat(3,i)],'g');
end

%% Calculate magnitude of N_hat
NN_hat = ppval(ppN_hat,xx);
NN_hat_mag = sqrt(sum((NN_hat.^2),1));

figure;
plot(xx,NN_hat_mag);
xlabel('x');
ylabel('Magnitude, N hat');