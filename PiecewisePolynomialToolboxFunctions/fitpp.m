function pp = fitpp(x,y,varargin)
% FITPP Fit a piecewise polynomial.
%   pp = FITPP(x,y) fits a piecewise polynomial to the x/y data provided.
%   The resultant polynomial is only C0 continuous without additional
%   conditions.
%
%   Data associated with a parametric equation can be specified by
%   adding rows to the provided y-data such that y(:,j) is taken an
%   m-dimensional array to be matched with x(j).
%       Example:
%           x = linspace(0,1,100);
%           y = rand(3,100);
%
%   pp = FITPP(x,y,x_d,dy) fits a piecewise polynomial to the x/y data
%   provided with constraints on the first derivative values dy specified 
%   at x_d.
%
%   pp = FITPP(x,y,x_d,dy,x_dd,ddy) fits a piecewise polynomial to the x/y
%   data provided with constraints on the first derivative values dy 
%   specified at x_d, and constraints on the second derivative values ddy 
%   specified at x_dd.
%
%   Derivative constraints can be added sequencially after the
%   initial data is provided. Constraints can be added to the 1st, 2nd,
%   3rd, 4th... derivative. Constraints on a given derivative can be
%   ignored using an empty set. For example:
%       pp = FITPP(x,y,[],[],x_dd,ddy) ignores the constraints on the first
%       derivative and applies provided contraints to the second
%       derivative.
%
%   --- NOTE: ADDITIONAL FIT PROPERTIES ARE NOT YET IMPLEMENTED ---
%   pp = FITPP(...,Name,Value) specifies additional fit properties.
%
%       Property Name | Property Values
%          Continuity | C0 - Ensures C0 continuity [Default]
%                     | C1 - Ensures C1 continuity
%                     | C2 - Ensures C2 continuity
%                     | ...
%                     | CN - Ensures CN continutiy
%   ---------------------------------------------------------------
%
%   See also spline, ppval, mkpp, unmkpp, diffpp, intpp
%
%   M. Kutzer, 17Nov2017, USNA

% Updates
%   02Feb2018 - Updated to correct for evolution of curve from 0 to arc
%               length
%   03Aug2018 - Updated to provide status updates for long processing

% TODO - Implement continuity constraints. 

%% Check inputs
% -> Check for minimum number of inputs
if nargin < 2
    error('Not enough input arguments.');
end
% -> Check for additional input pairs
if mod( numel(varargin),2 ) ~= 0
    error('Derivative constraints and properties must be specified in pairs.');
end
% -> Check that x-data is 1xN
dimX = size(x);
if numel(dimX) > 2 || dimX(1) > 1
    error('X-data must be specified by at most an 1xN array.');
end
% -> Check that y-data is MxN
% TODO - consider "two more values of than X has entries" note in spline.m
dimY = size(y);
if numel(dimY) > 2 || dimX(2) ~= dimY(2)
    error('Y-data must be specified by at most an MxN array given 1xN X-data.');
end

% TODO - Consider polyfun/private function chckxy.m

%% Order x and y
[x,idx] = sort(x,'ascend');
y = y(:,idx);

%% Parse inputs
% Set default continuity and initialize constraint variables
cN = 0;  % Differentiability class
dx = []; % x-values for derivative constraints
dy = []; % y-values for derivative constraints
dN = []; % derivative order

% TODO - this can be done much faster and in a cleaner way!
k = 0;
for i = 1:2:numel(varargin)
    if ischar( varargin{i} )
        % Assume property name is specified
        switch lower( varargin{i} )
            case 'continuity'
                % Define continuity
                cN = sscanf(lower(varargin{i+1}),'c%d',[1,1]);
                if isempty(cN)
                    error('fitPP:BadValue',...
                        'The property "%s" is not a correct value for "%s".',varargin{i+1},varargin{i});
                else
                    warning('Continuity constraints have not been implemented.');
                end
            otherwise
                error('fitPP:BadProperty',...
                    'The property "%s" is not recognized.',varargin{i});
        end
    else
        % Assume specific derivative constraints are provided
        k = k+1;
        if size(varargin{i},2) ~= size(varargin{i+1},2) || size(varargin{i+1},1) ~= dimY(1)
            error('Derivative constraints must be specified as a %dxK array given 1xK X-data.');
        end
        dx = [dx,varargin{i}];
        dy = [dy,varargin{i+1}];
        dN = [dN,repmat(k,1,numel(varargin{i}))];
    end
end

%% Fit initial set of polynomials

% Initialize coefficients array 
coefs = [];
% Initialize status flags
showWAITBAR = false;
showTIME = false;

nX = numel(x);
for i = 2:nX
    if i == 2
        ttt = tic;
    end
    
    % Find all user provided constraints
    bin = dx >= x(i-1) & dx <= x(i);
    const_x = [x((i-1):i),dx(bin)] - x(i-1); % "localize" x-values
    const_y = [y(:,(i-1):i),dy(:,bin)];      % y-values
    const_n = [0,0,dN(bin)];                 % derivative order
    
    % Define polynomial order
    % TODO - we should check if all constraints are unique
    n = numel(const_x)-1;
    
    % Define coefficients/powers
    P = zeros(n+1);
    for j = 0:n
        p = [(n-j):-1:0]';
        P(1:numel(p),j+1) = p;
    end
    C = [ones(size(P,1),1),P];
    
    % Define X
    X = [];
    % TODO - make this faster!
    for j = 1:numel(const_x)
        N = const_n(j);         % derivative order
        % Define powers
        if N > n
            p = zeros(n+1,1);
        else
            p = P(:,N+1);
        end
        % Define coefficients
        if N > n+1
            c = zeros(n+1,1);
        else
            c = prod(C(:,1:(N+1)),2);
        end
        % Calculate X
        % -> y = coef*X
        % X(:,j) = c.*(const_x(j).^p);
        tmpX = c.*(const_x(j).^p);
        if size(tmpX,1) > size(X,1)
            % -> New X is higher order, append zeros
            X_tmp = zeros(size(tmpX,1)-size(X,1),size(X,2));
            X = [X_tmp; X];
        elseif size(X,1) > size(tmpX,1)
            % -> Existing X is higher order, append zeros
            X_tmp = zeros(size(X,1)-size(tmpX,1),size(tmpX,2));
            tmpX = [X_tmp; tmpX];
        end
        X(:,j) = tmpX;
    end
    % Calculate coefficients
    % -> y = coef*X
    % -> coef = y*inv(X)
    coef = const_y/X; %y*inv(X);
    if size(coef,2) > size(coefs,2)
        % More new coefficients than existing coefficients
        % -> Pre-append zeros on existing coefficients
        coefs_tmp = zeros(size(coefs,1),size(coef,2)-size(coefs,2));
        coefs = [coefs_tmp,coefs];
        coefs = [coefs; coef];
    elseif size(coef,2) < size(coefs,2)
        % More existing coefficients than new coefficients
        coef_tmp = zeros(size(coef,1),size(coefs,2)-size(coef,2));
        % -> Pre-append zeros on new coefficients
        coef = [coef_tmp,coef];
        coefs = [coefs; coef];
    elseif size(coef,2) == size(coefs,2)
        % Same number of coefficients 
        coefs = [coefs; coef];
    end
    
    % --- Show fit status for long duration fits --------------------------
    if i == 2
        tSECONDS = toc(ttt)*nX;
        tMINUTES = tSECONDS/60;
        tHOURS   = tMINUTES/60;
        tDAYS    = tHOURS/24;
        
        if tSECONDS > 10
            showWAITBAR = true;
            hWAITBAR = waitbar(0,'Fitting piecewise polynomial...');
        end
        
        if tMINUTES > 0
            showTIME = true;
            fprintf('Status update for "%s.m"\n',mfilename);
            fprintf('               Current time: %s\n',datestr(now));
            fprintf('  Estimated completion time: %s\n',datestr(now + tDAYS));
        end
    end
    
    if showTIME
        if mod((i-1),4) == 0
            fprintf(char([8,8,8]));
        else
            fprintf('.');
        end
    end
    
    if showWAITBAR
        % Status update
        if ~ishandle(hWAITBAR)
            pp = [];
            fprintf( char(repmat(8,1,mod((i-1),4))) );
            warning('Fitting cancelled by user.');
            return
        end
        
        waitbar(i/nX,hWAITBAR);
    end
    % ---------------------------------------------------------------------
end

% --- Show fit status for long duration fits ------------------------------
if showWAITBAR
    delete(hWAITBAR);
end
if showTIME
    fprintf( char(repmat(8,1,mod((i-1),4))) );
    fprintf('     Actual completion time: %s\n',datestr(now));
end
% -------------------------------------------------------------------------

%% Apply CN constraint
% TODO - apply CN constraint

%% Package output
% Define breaks
breaks = x;
% Define dimension
d = dimY(1);
% Make piecewise polynomial
pp = mkpp(breaks,coefs,d);