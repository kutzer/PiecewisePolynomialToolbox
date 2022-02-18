function [tf,msg] = ispp(pp)
% ISPP returns true if the input is a valid piecewise polynomial structured
% array.
%   tf = ISPP(pp)
%
%   [tf,msg] = ISPP(___)
%
%   Input(s)
%       pp - MATLAB variable
%
%   Output(s)
%       tf  - binary scalar set to true if pp is a valid piecewise
%             polynomial structured array and false otherwise
%       msg - string argument describing array
%
%   See also fitpp intpp diffpp normalpp plotpp spline makepp unmkpp
%
%   M. Kutzer, 22Nov2017, USNA

% Updates
%   18Feb2022 - improved code speed

%% Check input(s)
narginchk(1,1);

%% Check variable
if ~isstruct(pp)
    % Set flag
    tf = false;
    % Create message
    if nargout > 1
        msg = sprintf('Specified variable is not a structured array');
    end

    return
end

% Specify required fields
fields = {...
    'form';...
    'breaks';...
    'coefs';...
    'pieces';...
    'order';...
    'dim'};
% Specify required form
form = 'pp';

% Check fields
bin = isfield(pp,fields);
if nnz(bin) == numel(fields)

    % Check for proper form
    formChk = {pp.form};
    bin = matches(formChk,form);
    if nnz(bin) ~= numel(pp)
        % Set flag
        tf = false;
        % Create message
        if nargout > 1
            badIdxs = find(~bin);
            msg = sprintf('The following elements of the specified variable are not form "%s": [',form);
            for i = 1:numel(badIdxs)
                msg = sprintf('%s%d',msg,badIdxs(i));
                if i ~= numel(badIdxs)
                    msg = sprintf('%s, ',msg);
                else
                    msg = sprintf('%s]',msg);
                end
            end
        end
        
        return
    end

    % Set flag
    tf = true;
    % Create message
    if nargout > 1
        if numel(pp) == 1
            msg = 'pp';
        else
            msg = 'ppArray';
        end
    end
else
    % Set flag
    tf = false;
    % Create message
    if nargout > 1
        missingFields = fields(~bin);
        msg = sprintf('Specified variable is missing the following fields: {');
        for i = 1:numel(missingFields)
            msg = sprintf('%s%s',msg,missingFields{i});
            if i ~= numel(missingFields)
                msg = sprintf('%s, ',msg);
            else
                msg = sprintf('%s}',msg);
            end
        end
    end
end