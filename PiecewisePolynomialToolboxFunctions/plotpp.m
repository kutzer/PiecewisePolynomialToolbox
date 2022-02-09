function [fig,axs] = plotpp(pp,t,varargin)
% PLOTPP plots an n-dimensional piecewise polynomial as a function of the
% parameterization term.
%   fig = plotpp(pp,t)
%
%   fig = plotpp(pp,t,n)
%
%   fig = plotpp(pp,t,n,lims)
%
%   [fig,axs] = plotpp(___)
%
%   Input(s)
%       pp   - piecewise polynomial structure
%       t    - array containing discrete points for evaluation
%       n    - [OPTIONAL] number of derivatives to evaluate in subplots
%       lims - [OPTIONAL] limits for polynomial (for display)
%           lims.pp(i,:)   - ith dimension limits of function, [lower, upper]
%           lims.dpp(i,:)  - ith dimension limits of first derivative, [lower, upper]
%           lims.ddpp(i,:) - ith dimension limits of second derivative, [lower, upper]
%           Etc.
%
%   Output(s)
%       fig - figure handle containing plots
%       axs - array of axes handles
%           axs(1) - original function
%           axs(2) - 1st derivative
%           axs(3) - 2nd derivative
%           Etc.
%
%   M. Kutzer, 08Feb2022, USNA

%% Parse input(s)
if nargin < 4
    lims = [];
else
    lims = varargin{2};
end

if nargin < 3
    n = 0;
else
    n = varargin{1};
end

% TODO - check inputs

%% Evaluate function(s)
di_pp = pp;
di_q{1} = ppval(di_pp,t);
for i = 1:n
    di_pp = diffpp(di_pp);
    di_q{i+1} = ppval(di_pp,t);
end

%% Parse breakpoints times and create xticklabels
t_b = pp.breaks;

xtck = t_b;
xtcklbl{1} = '$t_0$';
for i = 2:(numel(t_b)-1)
    xtcklbl{i} = sprintf('$t_{b_%d}$',i);
end
xtcklbl{i+1} = '$t_f$';

%% Define y-labels 
ylbl{1} = '$f(t)$';
ylbl{2} = '$\frac{df}{dt}(t)$';
for i = 2:n
    ylbl{i+1} = sprintf('$\\frac{d^{%d}f}{dt^{%d}}(t)$',i,i);
end

ytck{1} = [];
ytcklbl{1} = {};
for i = 1:n
    ytck{i+1} = 0;
    ytcklbl{i+1} = {'$0$'};
end

%% Create figure, axes, etc.
% Create figure
fig = figure('Color',[1 1 1],...
    'Units','Inches','Position',[0.25,0.65,6.27,5.40]);

for i = 1:(n+1)
    % Create axes
    axs(i) = subplot((n+1),1,i,'Parent',fig,'TickLabelInterpreter','Latex',...
        'XTickLabelMode','manual','XTickMode','manual',...
        'XTick',xtck,'XTickLabel',xtcklbl,...
        'YTickLabelMode','manual','YTickMode','manual',...
        'YTick',ytck{i},'YTickLabel',ytcklbl{i});
    %xlabel(axs(i),'$t$','Interpreter','Latex');
    ylabel(axs(i),ylbl{i},'Interpreter','Latex');
    hold(axs(i),'on');
    grid(axs(i),'on');
end

%% Define colors
colors = 'rgbcmyk';

%% Plot limits
% Initialize field names
fld{1} = 'pp';
for i = 1:n
    d = repmat('d',1,i);
    fld{i+1} = sprintf('%spp',d);
end

if ~isempty(lims)
    for i = 1:(n+1)
        if isfield(lims,fld{i})
            lim = lims.(fld{i});
            for j = 1:size(lim,1)
                w = 2 - 0.3*j;
                plim(i,j,1) = plot(axs(i),[t(1),t(end)],[lim(j,1),lim(j,1)],[colors(j),'--'],'LineWidth',w);
                plim(i,j,2) = plot(axs(i),[t(1),t(end)],[lim(j,2),lim(j,2)],[colors(j),'--'],'LineWidth',w);
                c = get(plim(i,j,1),'Color');
                c = c + 0.5*ones(1,3);
                c(c>1) = 1;
                set(plim(i,j,:),'color',c);
            end
        end
    end
end

%% Plot polynomial(s)
for i = 1:(n+1)
    for j = 1:size(di_q{i},1)
        w = 2 - 0.3*j;
        for k = 2:numel(t_b)
            bin = ( t >= t_b(k-1) & t < t_b(k) );
            t_i = t( bin );
            q_i = di_q{i}(j,bin);
            plt(i,j,k) = plot(axs(i),t_i,q_i,colors(j),'LineWidth',w);
        end
    end
end


%% Adjust axis limits
for j = 1:numel(axs)
    % -> Original axis limits
    xx0(1,:) = xlim(axs(j));
    yy0(1,:) = ylim(axs(j));
    % -> Buffered tight axis limits
    axis(axs(j),'tight'); drawnow;
    xx0 = xlim(axs(j));
    yy0 = ylim(axs(j));
    xx0 = (diff(xx0)/20)*[-1,1] + xx0;
    yy0 = (diff(yy0)/20)*[-1,1] + yy0;
    %xlim(axs(j),xx0);
    ylim(axs(j),yy0);
end