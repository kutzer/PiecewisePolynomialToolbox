function ipp = intpp(pp,c,x)
% INTPP Integrate piecewise polynomial
%   ipp = INTPP(pp) returns a piecewise polynomial associated with the
%   integral of pp assuming an initial condition of 0 at the first break.
%
%   ipp = INTPP(pp,c) returns a piecewise polynomial associated with the
%   integral of pp assuming an initial condition of "c" at the first break.
%
%   ipp = INTPP(pp,c,x) returns a piecewise polynomial associated with the
%   integral of pp assuming an initial condition of "c" when the integrated
%   function is evaluated at "x".

%   NOTE: This function calculates integration constants between breaks to
%   eliminate possible step discontinuities in the integral. The resultant
%   function should be C0 continuous.
%
%   M. Kutzer, 17Nov2017, USNA

%% Check number of input arguments
narginchk(1,3);

%% Extract breaks and coefficients from piecewise polynomial
[breaks,coeffs,~,~,d] = unmkpp(pp);

%% Set default values
if nargin < 2
    c = zeros(d,1);
end

if nargin < 3
    x = breaks(1);
end

%% Check initial condition
if numel(c) ~= d
    error('The number of initial conditions must match the dimension of the piecewise polynomial.');
end

%% Integrate piecewise polynomial (initial)
% Differentiate piecewise polynomial
%     f(x) =     a*x^(n)   +       b*x^(n-1) +       c*x^(n-2) ... 
% df(x)/dx = (n)*a*x^(n-1) + (n-1)*b*x^(n-2) + (n-2)*c*x^(n-3) ...
icoeffs = coeffs * diag( [size(coeffs,2):-1:1] )^(-1);
icoeffs(:,end+1) = 0;

%% Impose C0
% Define polynomial order
n = size(icoeffs,2) - 1;

idx0 = 1:d:size(icoeffs,1);
idxF = d:d:size(icoeffs,1);
for i = 2:numel(breaks)-1
    XX = (breaks(i) - breaks(i-1)).^[n:-1:0]';
    cc = icoeffs(idx0(i-1):idxF(i-1),:)*XX;
    icoeffs(idx0(i):idxF(i),end) = cc;
end

%% Impose "c" and "x"
ipp = mkpp(breaks,icoeffs,d);
cc = ppval(ipp,x);
c_offset = c - cc;

icoeffs(:,end) = icoeffs(:,end) + repmat(c_offset,size(icoeffs,1)/d,1);

%% Create final integrated piecewise polynomial
ipp = mkpp(breaks,icoeffs,d);
