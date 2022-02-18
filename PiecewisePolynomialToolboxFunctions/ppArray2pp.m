function pp = ppArray2pp(ppArray)
% PPARRAY2PP converts an n-element structured array containing piecewise
% polynomial definitions into a single piecewise polynomial.
%   pp = ppArray2pp(ppArray)
%
%   Input(s)
%       ppArray - n-element piecewise polynomial structured array
%
%   Output(s)
%       pp - m-dimensional piecewise polynomial structured array
%
%   NOTE: This function stacks the piecewise polynomial dimensions in order
%         and, if breakpoints are not shared, creates a common set of
%         shared breakpoints for the resultant function.
%
%   Example:
%     %% Define points & independent variable for 1D piecewise fit
%     x1 = rand(1,5);                   % 1D points
%     s1 = linspace(0,1,size(x1,2));    % Independent variable
%     %% Fit piecewise polynomial (cubic spline)
%     ppArray(1) = spline(s1,x1);
%     %% Define points & independent variable for 3D piecewise fit
%     x2 = rand(3,4);                   % 3D points
%     s2 = linspace(-1,2,size(x2,2));   % Independent variable
%     %% Fit piecewise polynomial (C^0 continuous)
%     ppArray(2) = fitpp(s2,x2);
%     %% Combine piecewise polynomial array to a single piecewise polynomial
%     pp = ppArray2pp(ppArray);
%
%   M. Kutzer, 18Feb2022, USNA

pp = [];

%% Check input(s)
narginchk(1,1);

[tf,msg] = ispp(ppArray);
if ~tf
    error('Input variable must be a valid piecewise polynomial structured array\n\n%s',msg);
end

%% Check for single element ppArray
if numel(ppArray) == 1
    % Input piecewise polynomial array is a single piecewise polynomial
    pp = ppArray;
    return
end

%% Unmake each piecewise polynomial array element
for i = 1:numel(ppArray)
    [breaks_i{i},coefs_i{i},L_i(i),order_i(i),dim_i(i)] = unmkpp(ppArray(i));
end

%% Combine breakpoints
breaks = sort(unique([breaks_i{:}]));

% Check initial and final break points
% NOTE: The following syntax is not supported
%   bin0 = breaks(1) == [breaks_i{:}(1)]
%   binf = breaks(end) == [breaks_i{:}(end)]
for i = 1:numel(breaks_i)
    breaks_i_0(i) = breaks_i{i}(1);
    breaks_i_f(i) = breaks_i{i}(end);
end
bin_0 = breaks(1) == breaks_i_0;
bin_f = breaks(end) == breaks_i_f;

% Check initial breakpoint
if nnz(bin_0) ~= numel(bin_0)
    warning([...
        'The specified piecewise polynomials do not share the same initial breakpoint.\n',...
        'The initial breakpoint will be changed from %.4f to the first overlapping \n',...
        'breakpoint of %.4f.'],breaks(1),max(breaks_i_0));
    bin = breaks >= max(breaks_i_0);
    breaks = breaks(bin);
end

% Check final breakpoint
if nnz(bin_f) ~= numel(bin_f)
    warning([...
        'The specified piecewise polynomials do not share the same final breakpoint.\n',...
        'The final breakpoint will be changed from %.4f to the final overlapping \n',...
        'breakpoint of %.4f.'],breaks(end),min(breaks_i_f));
    bin = breaks <= min(breaks_i_f);
    breaks = breaks(bin);
end

%% Make coefficients a common order
order = max(order_i);
for i = 1:numel(coefs_i)
    m = size(coefs_i{i},1);
    coefs_c{i} = zeros(m,order);
    coefs_c{i}(:,(end-order_i(i)+1):end) = coefs_i{i};
end
coefs_i = coefs_c;
clear coefs_c

%% Isolate coefficients by dimensions
% Pattern for native coefficients:
% coefs_i{i}:
%   [dim_1, breakpoint_i_1 to breakpoint_i_2]
%   [dim_2, breakpoint_i_1 to breakpoint_i_2]
%   [ ... , breakpoint_i_1 to breakpoint_i_2]
%   [dim_d, breakpoint_i_1 to breakpoint_i_2]
%   [dim_1, breakpoint_i_2 to breakpoint_i_3]
%   [dim_2, breakpoint_i_2 to breakpoint_i_3]
%   [ ... , breakpoint_i_2 to breakpoint_i_3]
%   [dim_d, breakpoint_i_2 to breakpoint_i_3]
%   [ ... ,  ...           to  ...          ]
%   [dim_1, breakpoint_i_{f-1} to breakpoint_i_f]
%   [dim_2, breakpoint_i_{f-1} to breakpoint_i_f]
%   [ ... , breakpoint_i_{f-1} to breakpoint_i_f]
%   [dim_d, breakpoint_i_{f-1} to breakpoint_i_f]
%
% Common pattern
% coefs_c{i}{j}
%   [dim_j, breakpoint_i_1 to breakpoint_i_2]
%   [dim_j, breakpoint_i_2 to breakpoint_i_3]
%   [ ... ,  ...           to  ...          ]
%   [dim_j, breakpoint_i_{f-1} to breakpoint_i_f]
for i = 1:numel(coefs_i)
    n = size(coefs_i{i},1);
    for d = 1:dim_i(i)
        idx = d:dim_i(i):n;
        coefs_c{i}{d} = coefs_i{i}(idx,:);
    end
end

%% Create common set of coefficients for new breakpoints
c(:,order) = ones(order,1);
for n = order:-1:1
    p(:,n) = zeros(order,1);
    p(1:n,n) = [(n-1):-1:0].';

    if n > 1
        c(:,n-1) = c(:,n).*p(:,n);
    end
end
%p = p + tril(ones(size(p)),-1);
vs = @(s,i)c(:,i).*( s.*ones(order,1) ).^( p(:,i) );
vS = @(s)c.*( s.*ones(order,order) ).^( p );
for i = 1:numel(coefs_c)
    for d = 1:numel(coefs_c{i})
        % Initialize common breakpoints coefficients
        coefs_n{i}{d} = [];
        for b = 2:numel(breaks_i{i})
            % Extract current polynomial coefficients
            C = coefs_c{i}{d}(b-1,:);

            % Isolate breakpoints for current polynomial
            s0_i = breaks_i{i}(b-1); % Current initial breakpoint
            s1_i = breaks_i{i}(b);   % Current final breakpoint
            
            % Check for overlaps with new breakpoints
            bin0_i = s0_i == breaks;
            bin1_i = s1_i == breaks;

            % WE HAVE NOT ADDRESSED THE CASE WHERE A FUNCTION HAS BROADER
            % BREAKPOINTS THAN THE COMBINE BREAKPOINTS! 
            binTEST = s0_i <= breaks & s1_i >= breaks % <--- USE THIS DUMMY!

            bin = s0_i == breaks | s1_i == breaks;

            if nnz(bin) == 0
                % Breakpoint interval is not included in combine pp
                fprintf('IGNORE\n')
                fprintf('%d, %d [%.4f,%.4f]: %d\n',i,d,s0_i,s1_i,nnz(bin));
                disp(bin)
                continue
            end

            if nnz(bin) == 1
                % Special case?!?!
            end

            % Differentiate overlap condition
            idx = find(bin);
            nParts = diff(idx);

            if nParts == 1
                fprintf('KEEP AS-IS\n')
                fprintf('%d, %d [%.4f,%.4f]: %d\n',i,d,s0_i,s1_i,nnz(bin));
                disp(bin)
                % Append coefficients
                coefs_n{i}{d}(end+1,:) = C;
                continue
            end

            if nParts > 1
                fprintf('SPLIT INTO %d parts\n',nParts)
                fprintf('%d, %d [%.4f,%.4f]: %d\n',i,d,s0_i,s1_i,nnz(bin));
                disp(bin)

                % Define initialize condition for shared initial breakpoint
                for b_n = (idx(1)+1):idx(2)
                    s0_n = breaks(b_n-1); % Initial breakpoint
                    s1_n = breaks(b_n);   % Final breakpoint

                    % Update initial conditions for upcoming initial breakpoint
                    % -> Preserve C^(order-1) continuity
                    C0_n = C*vS(s0_n - s0_i);
                    %C1_n = C*vS(s1_n - s0_i); % <--- Unused
                    
                    % Append coefficients
                    C_n = C0_n;
                    coefs_n{i}{d}(end+1,:) = C_n;

                    % DEBUG
                    % -> Evaulate existing coeffients
                    f0_i = C*vS(s0_n - s0_i);
                    f1_i = C*vS(s1_n - s0_i);
                    % -> Evaulate new coefficients
                    f0 = C_n*vS(s0_n - s0_n);
                    f1 = C_n*vS(s1_n - s0_n);

                    % DEBUG
                    ff0_i = fliplr(f0_i);
                    ff1_i = fliplr(f1_i);
                    fprintf('s   = {%9.4f,%9.4f}, ',s0_n,s1_n);
                    for k = 1:order
                        strds = repmat('d',1,k-1);
                        fprintf('%sf   = {%9.4f,%9.4f}, ',strds,ff0_i(k),ff1_i(k))
                    end
                    fprintf('\n');

                    ff0 = fliplr(f0);
                    ff1 = fliplr(f1);
                    fprintf('s_%d = {%9.4f,%9.4f}, ',b_n-1,s0_n,s1_n);
                    for k = 1:order
                        strds = repmat('d',1,k-1);
                        fprintf('%sf_%d = {%9.4f,%9.4f}, ',strds,b_n-1,ff0(k),ff1(k));
                    end
                    fprintf('\n');
                end
                fprintf('\n')
            end
        end
    end
end

%% Combine ordered coefficients
L = numel(breaks)-1;    % Number of intervals
dim = sum(dim_i);
coefs = [];
for k = 1:L
    for i = 1:numel(coefs_n)
        for d_i = 1:numel(coefs_n{i})
            coefs(end+1,:) = coefs_n{i}{d_i}(k,:);
        end
    end
end

%% Make combine piecewise polynomial
pp = mkpp(breaks,coefs,dim);