function varargout = normalpp(varargin)
% NORMALPP Calculate a function for the normal to a piecewise polynomial.
%   N = NORMALPP(pp) returns a piecewise polynomial associated with the
%   normal of pp with a default of 50 spline segments between break points
%   when defining N.
%
%   N = NORMALPP(pp,n) returns a piecewise polynomial associated with the
%   normal of pp with n spline segments between each break point when 
%   defining N.
%
%   [N, N_hat] = NORMALPP(pp) returns piecewise polynomials associated
%   with the normal and unit normal of pp with a default of 50 spline
%   segments between break points when defining N, and an additional 0  
%   spline segments between break points of N when defining N-hat.
%
%   [N, N_hat] = NORMALPP(pp,n) returns piecewise polynomials associated
%   with the normal and unit normal of pp with n spline segments between 
%   break points when defining N, and an additional 0 spline segments 
%   between break points of N when defining N-hat.
%
%   [N, N_hat] = NORMALPP(pp,n,m) returns piecewise polynomials associated
%   with the normal and unit normal of pp with n spline segments between 
%   break points when defining N, and an additional m spline segments 
%   between break points of N when defining N-hat.
%
%   See also spline ppval mkpp unmkpp fitpp diffpp intpp plotpp tangentpp
%   framepp arcLengthParamPP appendpp ispp ppArray2pp
%
%   M. Kutzer, 17Jan2018, USNA

% Updates
%   20Feb2022 - Added "See also"

%% Check number of input arguments
narginchk(1,3);

pp = varargin{1};
if nargin < 2
    n = 50;
else
    n = varargin{2};
end

if nargin < 3
    m = 0;
else
    m = varargin{3};
end

%% Define unit tangent
[~,T_hat] = tangentpp(pp,n);

%% Differentiate piecewise polynomial of unit tangent to define normal
N = diffpp(T_hat);

%% Package output
varargout{1} = N;
if nargout == 1
    return
end

%% Calculate normal
breaks = unmkpp(N);
x = breaks;

xx = [];
for i = 2:numel(x)
    xx_tmp = linspace(x(i-1),x(i),m+2);
    xx_tmp(end) = [];
    xx = [xx, xx_tmp];
end
xx(end+1) = x(end);

Y = ppval(N,xx);

Y_hat = Y ./ sqrt(sum((Y.^2),1));

N_hat = spline(xx,Y_hat);
varargout{2} = N_hat;
