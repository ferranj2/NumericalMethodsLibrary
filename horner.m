%% Horner's multiplication rule
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE), College of Arts and Sciences (COAS)
%  For use in MA 412, MA 413, AE 435, AE 440 and any other course that
%  would benefit from a data-fitting tool.
%% Description
% Evaluates a polynomial specified by an array of coefficients at specified
% values of the indeterminate variable. The size of the array must be of
% size "n+1", where "n" is the order of the polynomial. Horner's rule
% evaluates a polynomial using "n" additions and multiplications.
%% Formulae
% $P_{n}(x) = a_{0} + a_{1}x + a_{2}x^{2} + ... + a_{n}x^{n}$
%%
% $P_{n}(x) = a_{0} + x(a_{1} + x(a_{2} + x(a_{3}+)...x(a_{n-1} + a_{n}x)...)))$
%% Required Plugins
% * none
%% Changelog
%  v1.0,(08/07/2022): Initial Release. Reinvented the wheel once again.
%% Syntax
% * INPUT(*coeff*): A row or column array that lists the coefficients of
% the powers of the indeterminates in ascending order (i.e., the constant
% coefficient comes first, followed by the coefficients of the ascending
% powers).
% * INPUT(*x*): Array of numeric values of the indeterminate to evaluate
% polynomial at.
% * OUTPUT(*val*): Values of the evaluated polynomial at the points
% specified by *x*.
%% Function definition
function val = horner(coeff,x)
[rc,cc] = size(coeff);
if rc ~= 1 && cc ~= 1 %
    error('Must input a column/row array for the coefficients!');
end
n_coeff = length(coeff);
val = ones(size(x))*coeff(n_coeff); %Vector initialized to last coefficient.
if n_coeff == 1 %Input is a scalar.
    return;
end
val = val.*x;
for ii = n_coeff-1:-1:2
    val = (val + coeff(ii)).*x;
end
val = val + coeff(1);
end