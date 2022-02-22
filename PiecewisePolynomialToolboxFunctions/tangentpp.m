function varargout = tangentpp(varargin)
% TANGENTPP Calculate a function for the tangent of a piecewise
% polynomial.
%   T = TANGENTPP(pp) returns a piecewise polynomial associated with the
%   tangent of pp.
%
%   [T, T_hat] = TANGENTPP(pp) returns piecewise polynomials associated
%   with the tangent and unit tangent of pp with a default of 50 spline 
%   segments between break points when defining T_hat. 
%
%   [T, T_hat] = TANGENTPP(pp,N) returns piecewise polynomials associated
%   with the tangent and unit tangent of pp with a default of N spline 
%   segments between break points when defining T_hat. 
%
%   See also spline ppval mkpp unmkpp fitpp diffpp intpp plotpp 
%   normalpp framepp arcLengthParamPP appendpp ispp ppArray2pp
%
%   M. Kutzer, 17Jan2018, USNA

% Updates
%   20Feb2022 - Added "See also"

%% Check number of input arguments
narginchk(1,2);

pp = varargin{1};
if nargin < 2
    N = 50;
else
    N = varargin{2};
end

%% Differentiate piecewise polynomial to define tangent
T = diffpp(pp);

%% Package output
varargout{1} = T;
if nargout == 1
    return
end

%% Calculate unit tangent
breaks = unmkpp(T);
x = breaks;

xx = [];
for i = 2:numel(x)
    xx_tmp = linspace(x(i-1),x(i),N+2);
    xx_tmp(end) = [];
    xx = [xx, xx_tmp];
end
xx(end+1) = x(end);

Y = ppval(T,xx);

Y_hat = Y ./ sqrt(sum((Y.^2),1));

T_hat = spline(xx,Y_hat);
varargout{2} = T_hat;