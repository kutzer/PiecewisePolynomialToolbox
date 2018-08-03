function tf = ispp(pp)
% ISPP True for valid piecewise polynomial structured arrays.
%   tf = ISPP(pp) returns an array that contains 1's where the elements of
%   pp are valid piecewise polynomial structured arrays.
%
%   See also makepp unmkpp
%
%   M. Kutzer, 22Nov2017, USNA

%% Check number of inputs
narginchk(1,1);

%% Define field names 
tf = false;
fnames = {'form'; 'breaks'; 'coefs'; 'pieces'; 'order'; 'dim'};
if isstruct(pp)
    % TODO - consider checking all fields
    if sum( isfield(pp,fnames) ) == numel(fnames)
        if strcmp(pp.form,'pp')
            tf = true;
        end
    end
end