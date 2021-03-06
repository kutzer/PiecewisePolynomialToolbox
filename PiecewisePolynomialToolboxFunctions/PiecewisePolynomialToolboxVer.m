function varargout = PiecewisePolynomialToolboxVer
% PIECEWISEPOLYNOMIALTOOLBOXVER displays the Piecewise Polynomial Toolbox information.
%   PIECEWISEPOLYNOMIALTOOLBOXVER displays the information to the command prompt.
%
%   A = PIECEWISEPOLYNOMIALTOOLBOXVER returns in A the sorted struct array of  
%   version information for the Piecewise Polynomial Toolbox.
%     The definition of struct A is:
%             A.Name      : toolbox name
%             A.Version   : toolbox version number
%             A.Release   : toolbox release string
%             A.Date      : toolbox release date
%
%   M. Kutzer 05Feb2018, USNA

A.Name = 'Piecewise Polynomial Toolbox';
A.Version = '1.0.3';
A.Release = '(R2019b)';
A.Date = '08-Jan-2021';
A.URLVer = 1;

msg{1} = sprintf('MATLAB %s Version: %s %s',A.Name, A.Version, A.Release);
msg{2} = sprintf('Release Date: %s',A.Date);

n = 0;
for i = 1:numel(msg)
    n = max( [n,numel(msg{i})] );
end

fprintf('%s\n',repmat('-',1,n));
for i = 1:numel(msg)
    fprintf('%s\n',msg{i});
end
fprintf('%s\n',repmat('-',1,n));

if nargout == 1
    varargout{1} = A;
end