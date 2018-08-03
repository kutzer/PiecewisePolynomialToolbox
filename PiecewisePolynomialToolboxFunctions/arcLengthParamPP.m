function varargout = arcLengthParamPP(varargin)
% ARCLENGTHPARAMPP Arc-length parameterization of piecewise polynomial(s).
%   pps = arcLengthParamPP(pp) returns the arc length parameterized curve 
%   with a default of 50 spline segments between break points.
%
%   pps = arcLengthParamPP(pp,N) returns the arc length parameterized curve 
%   with a user defined N spline segments between break points.
%
%   [pps,slim] = arcLengthParamPP(___) returns the arc length parametrized
%   curve and the limits of arc length [0, s_max].
%
%   See also spline fitpp
%
%   M. Kutzer, 13Dec2017, USNA

%% Check number of inputs
narginchk(1,2);

%% Define defaults
pp = varargin{1};
if ~ispp(pp)
    error('Input must be a valid piecewise polynomial structure.');
end

if nargin < 2
    N = 50;
else
    N = varargin{2};
end
    
%% Fit arc length s as a function input x
% -> Define tangent
% pp  - f(x)
% ppT - T(x) = df(x)/dx
ppT = diffpp(pp);

% -> Extract breaks from piecewise polynomial
[breaks] = unmkpp(pp);

% Define new breaks
n = numel(breaks);
newBreaks = zeros(1,N*(n-1)+1);
for i = 2:numel(breaks)
    % Define N intermediate break points
    tmpBreaks = linspace(breaks(i-1),breaks(i),N+1);
    tmpBreaks(end) = [];
    % Combine to new break points
    newBreaks(1,(N*(i-2))+(1:N)) = tmpBreaks;
end
% Add final break point
newBreaks(end) = breaks(end);

% -> Evaluate tangent at new breakpoints
T = ppval(ppT,newBreaks);
% -> Find and fit norm of tangent
normT = sqrt( sum(T.^2,1) );
% TODO - Consider analytically finding dNormT and using fitpp.m
% ppNormT - ds(x)/dx = |T(x)|
ppNormT = spline(newBreaks,normT);
% -> Integrate norm of tangent to find and evaluate s(x)
c0 = 0;              % Initial condition at first break point
x0 = newBreaks(1);   % First breakpoint
% ppS_X - s(x)
ppS_X = intpp(ppNormT,c0,x0);
s_x = ppval(ppS_X,newBreaks);

%% Find function input x as a function of arc length s
% -> Estimate and fit inverse of s(x) to find x(s)
% ppX_S - x(s), inverse of s(x)
ppX_S = spline(s_x,newBreaks);

% -> Uniformly redistribute s(x) to evaluate x(s)
% TODO - See [1] for a better way of looking at this fit:
%   [1] H. Wang, J. Kearney, K. Atkinson, "Arc-Length Parameterized Spline
%   Curves for Real-Time Simulation"
%
% --> Current method feels like resampling to hand-wave our way out of a
%     redundant step.
s = linspace(min(s_x),max(s_x),numel(s_x));
x_s = ppval(ppX_S,s); 

%% Create arc-length parameterization
pps = fitpp(s,ppval(pp,x_s),s,ppval(ppT,x_s)./ppval(ppNormT,x_s));

%% Package outputs
if nargout > 0
    varargout{1} = pps;
end
if nargout > 1
    varargout{2} = [min(s),max(s)];
end