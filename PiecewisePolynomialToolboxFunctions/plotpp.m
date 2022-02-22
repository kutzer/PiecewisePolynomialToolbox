function hndls_s = plotpp(pp,t,varargin)
% PLOTPP plots an n-dimensional piecewise polynomial as a function of the
% parameterization term.
%   hndls_s = plotpp(pp,t)
%
%   hndls_s = plotpp(pp,t,n)
%
%   hndls_s = plotpp(pp,t,n,lims,indStr)
%
%   hndls_s = plotpp(pp,t,n,lims,indStr)
%
%   hndls_s = plotpp(pp,t,n,lims,indStr,fcnStr)
%
%   [hndls_s,hndls_f] = plotpp(___)
%
%   Input(s)
%       pp     - piecewise polynomial structure
%       t      - array containing discrete points for evaluation
%       n      - [OPTIONAL] number of derivatives to evaluate in subplots
%       lims   - [OPTIONAL] limits for polynomial (for display). Define
%                lims = [] to skip this option. 
%           lims.pp(i,:)   - ith dimension limits of function, [lower, upper]
%           lims.dpp(i,:)  - ith dimension limits of first derivative, [lower, upper]
%           lims.ddpp(i,:) - ith dimension limits of second derivative, [lower, upper]
%           Etc.
%       indStr - [OPTIONAL] string defining independent variable (default is 's') 
%       fcnStr - [OPTIONAL] string defining function variable (default is 'f')
%
%    fcnStr NOTE: For adding vector symbols or other latex '\*' symbols:
%                 Use '\\\' instead of '\'
%        EXAMPLE: fcnStr = '\vec{q}' should be fcnStr = '\\\vec{q}'
%
%   Output(s)
%       fig - figure handle containing plots
%       axs - array of axes handles
%           axs(1) - original function
%           axs(2) - 1st derivative
%           axs(3) - 2nd derivative
%           Etc.
%
% Example:
%   %% Define points & independent variable for 1D piecewise fit
%   x = [1,5,1,-8,2]; % Points
%   s = 1:numel(x);   % Independent variable
%   %% Fit piecewise polynomial (C^0 continuity)
%   pp = fitpp(s,x);
%   %% Plot piecewise polynomial with 1st & 2nd derivatives
%   ss = linspace(s(1),s(end),1000);
%   fig0 = plotpp(pp,ss,2,[],'s','f');
%
%   See also spline ppval mkpp unmkpp fitpp diffpp intpp tangentpp
%   normalpp framepp arcLengthParamPP appendpp ispp ppArray2pp
%
%   M. Kutzer, 08Feb2022, USNA

% Updates
%   18Feb2022 - Added fcnStr
%   20Feb2022 - Added "See also"

%% Parse input(s)
if nargin < 6
    fcnStr = 'f';
else
    fcnStr = varargin{4};
end

if nargin < 5
    indStr = 's';
else
    indStr = varargin{3};
end

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

%% Check input(s)
[tf,msg] = ispp(pp);
if ~tf
    error('Input variable "pp" must be a valid piecewise polynomial structured array\n\n%s',msg);
end

switch msg
    case 'pp'
        % Do nothing
    case 'ppArray'
        pp = ppArray2pp(pp);
    otherwise
        error('Unexpected "ispp" message.');
end

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
for i = 0:numel(t_b)
    xtcklbl{i+1} = sprintf('$%s_{%d}$',indStr,i);
end

%% Define y-labels 
ylbl{1} = sprintf('$%s(%s)$',fcnStr,indStr);
ylbl{2} = sprintf('$\\frac{d%s}{d%s}(%s)$',fcnStr,indStr,indStr);
for i = 2:n
    ylbl{i+1} = sprintf('$\\frac{d^{%d}%s}{d%s^{%d}}(%s)$',i,fcnStr,indStr,i,indStr);
end

ytck{1} = [];
ytcklbl{1} = {};
for i = 1:n
    ytck{i+1} = 0;
    ytcklbl{i+1} = {'$0$'};
end

%% ------------------------------------------------------------------------
%  ---- PLOT DEPENDENT VARIABLES vs INDEPENDENT VARIABLE ------------------
%  ------------------------------------------------------------------------
%% Create figure, axes, etc.
% Create figure
fig = figure('Color',[1 1 1],...
    'Units','Inches','Position',[0.25,0.65,6.27,5.40],'Name','plotpp: Dependend vs Independent');

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

% Initialize handle array
plim = nan(n+1,size(di_q{i},1),2);
size(plim)
% Populate handle array
if ~isempty(lims)
    for i = 1:(n+1)
        if isfield(lims,fld{i})
            lim = lims.(fld{i});
            for j = 1:size(lim,1)
                % Adjust width so all lines are visuble
                w = 2 - 0.3*j;
                % Plot limits
                plim(i,j,1) = plot(axs(i),[t(1),t(end)],[lim(j,1),lim(j,1)],[colors(j),'--'],'LineWidth',w);
                plim(i,j,2) = plot(axs(i),[t(1),t(end)],[lim(j,2),lim(j,2)],[colors(j),'--'],'LineWidth',w);
                % Make colors a shade lighter
                c = get(plim(i,j,1),'Color');
                c = c + 0.5*ones(1,3);
                c(c>1) = 1;
                % Update colors
                set(plim(i,j,:),'color',c);
            end
        end
    end
end
size(plim)

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
    xx0 = (diff(xx0)/100)*[-1,1] + xx0;
    yy0 = (diff(yy0)/20)*[-1,1] + yy0;
    xlim(axs(j),xx0);
    ylim(axs(j),yy0);
end

%% Package output
hndls_s.Figure = fig;
hndls_s.Axes = axs;
hndls_s.Plots = plt;
hndls_s.Limits = plim;

if nargout <= 1
    return
end

%% ------------------------------------------------------------------------
%  ---- PLOT DEPENDENT VARIABLES vs DEPENDENT VARIABLES -------------------
%  ------------------------------------------------------------------------

