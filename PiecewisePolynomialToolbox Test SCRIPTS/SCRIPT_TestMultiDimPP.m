%% SCRIPT_TestMultiDimPP
x = linspace(-2,2,5);

f{1}  = @(x_in)  2*x.^2 + 1;
f{2}  = @(x_in)  4*x.^1 + 2;
f{3}  = @(x_in) -6*x.^2 + 3;

for i = 1:numel(f)
    y(i,:) = f{i}(x);
    pp{i} = spline(x,y(i,:));
end

pp_all = spline(x,y);

for i = 1:3
    [breaks{i},coeffs{i},~,~,d] = unmkpp(pp{i});
    breaks{i}
    coeffs{i}
    d
end


[breaks_all,coeffs_all,~,~,d] = unmkpp(pp_all)