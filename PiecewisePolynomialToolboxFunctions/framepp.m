function F = framepp(varargin)
% FRAMEPP Define a piecewise polynomial function defining the position and 
% orientation of a frame along a prescribed piecewise polynomial.
%   F = framepp(pp) returns a piecewise polynomial function defining the
%   position and orientation of a frame along the prescribed function such
%   that
%       q = ppval(F,x) returns a 6-element vector containing position
%       q(1:3) and rotation expm(wedge(q(4:6)) information. 
%
%       -> The x-direction of the rotation is defined by the tangent
%       -> The y-direction of the rotation is defined by the normal
%
%   F = framepp(pp,N) specifies N added break points when defining the
%   unit tangent (see tangentpp.m).
%
%   F = framepp(pp,N,M) specifies N added break points when defining the
%   unit tangent (see tangentpp.m) and M added break points when defining
%   the unit normal (see normalpp.m).
%
%   M. Kutzer, 17Jan2018, USNA

%% Check number of input arguments
narginchk(1,3);

pp = varargin{1};
if nargin < 2
    N = 50;
else
    N = varargin{2};
end

if nargin < 3
    M = 0;
else
    M = varargin{3};
end

%% Define unit tangent
[~,T_hat] = tangentpp(pp,N);

%% Define unit normal
[~,N_hat] = normalpp(pp,N,M);

%% Define q to fit polynomial
breaks = unmkpp(N_hat);
x = breaks;

% Position - q(1:3,:)
q = ppval(pp,x);
[m,n] = size(q)

% Orientation
x_hat = ppval(T_hat,x);
y_hat = ppval(N_hat,x);

theta_all = [];
k_all = [];
for i = 1:numel(x)
    z_hat(:,i) = cross(x_hat(:,i),y_hat(:,i));
    R = [x_hat(:,i),y_hat(:,i),z_hat(:,i)];
    
    [k,theta] = SOtoAxisAngle(R);
    if i > 1
        out = dot(k_all(:,i-1),k);
    else
        out = 1;
    end
    
    if out > 0
        theta_all(end+1) = wrapTo2Pi(theta);
        k_all(:,end+1) = k;
    else
        theta_all(end+1) = wrapTo2Pi(-theta);
        k_all(:,end+1) = -k;
    end
    
    q((m+1):(m+numel(k)),i) = k_all(:,end)*theta_all(end);
end

%% Define piecewise polynomial
F = spline(x,q);
