function pp = appendpp(pp1,pp2,varargin)
% APPENDPP Appends two piecewise polynomials.
%   pp = APPENDPP(pp1,pp2) appends two piecewise polynomials by appending
%   the associated fields to create a combine piecewise polynomial. Unless
%   otherwise specified, a break offset between the two piecewise
%   polynomials equal to the mean difference of all breaks is used.
%
%   pp = APPENDPP(pp1,pp2,breakOffset) allows the user to specify the break
%   offset used to combine the two sets of break points.
%
%   See also spline ppval mkpp unmkpp fitpp diffpp intpp plotpp tangentpp
%   normalpp framepp arcLengthParamPP appendpp ispp ppArray2pp
%
%   M. Kutzer, 15Oct2019, USNA

% Updates
%   20Feb2022 - Added "See also"

%% Check number of input arguments
narginchk(2,3);

%% Unmake and combine piecewise polynomials
% Extract breaks, coefficients, number of intervals, order, and dimension
% from piecewise polynomials
[breaks{1},coeffs{1},L(1),k(1),dim(1)] = unmkpp(pp1);
[breaks{2},coeffs{2},L(2),k(2),dim(2)] = unmkpp(pp2);

% Define input variable names for warnings & errors
varNames = {'pp1','pp2'};

% Check dimension
if dim(1) ~= dim(2)
    error('The dimension of "%s" is %d and the dimension of "%s" is %d. These values must be the same to append these piecewise polynimials.',...
        varNames{1},dim(1),varNames{2},dim(2));
else
    dim = dim(1);   % One dimension to rule them all...
end
        
% Check degree 
% Match the order of the piecewise polynomials
if k(1) ~= k(2)
    
    [kMax,idxMax] = max(k);
    [kMin,idxMin] = min(k);
    
    % Warn user 
    warning('The order of "%s" is %d and the order of "%s" is %d. Updating "%s" to an order of %d.',...
        varNames{1},k(1),varNames{2},k(2),varNames{idxMin},kMax);
    
    % Update coefficients to match order
    coeffsTMP(:,(kMax-kMin+1):kMax) = coeffs{idxMin};
    coeffs{idxMin} = coeffsTMP;
    
    % Update order
    k(idxMin) = kMax;
end

%% Define break offset
if nargin < 3
    dbreaks = [];
    for i = 1:numel(breaks)
        dbreaks = [dbreaks, diff(breaks{i})];
    end
    breakOffset = mean(dbreaks);
else
    breakOffset = varargin{1};
end

%% Append breaks and coefficients
breaks = [breaks{1}, breaks{2} + repmat(breakOffset + breaks{1}(end), size(breaks{2}))];
coeffs = [coeffs{1}; coeffs{2}];

%% Make pp
pp = mkpp(breaks,coeffs,dim);