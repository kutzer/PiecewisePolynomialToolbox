function dpp = diffpp(pp)
% DIFFPP Differentiate a piecewise polynomial.
%   dpp = DIFFPP(pp) returns a piecewise polynomial associated with the
%   derivative of pp.
%
%   M. Kutzer, 17Nov2017, USNA

% Updates:
%   09Feb2022 - If the order of the polynomial fit falls below 1, the
%               coefficients of 0 are added to the first order term

%% Check number of input arguments
narginchk(1,1);

%% Differentiate piecwise polynomial
% Extract breaks and coefficients from piecewise polynomial
[breaks,coeffs,~,~,d] = unmkpp(pp);

% Differentiate piecewise polynomial
%     f(x) =     a*x^(n)   +       b*x^(n-1) +       c*x^(n-2) ...
% df(x)/dx = (n)*a*x^(n-1) + (n-1)*b*x^(n-2) + (n-2)*c*x^(n-3) ...
dcoeffs = coeffs * diag( [(size(coeffs,2)-1):-1:0] );
dcoeffs(:,end) = [];

% Maintain at least a first order polynomial
if size(dcoeffs,2) < 2
    dcoeffs = [zeros(size(dcoeffs,1)), dcoeffs];
end

%% Create the piecewise polynomial
dpp = mkpp(breaks,dcoeffs,d);