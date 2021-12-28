%% Polynomial Numeric evaluation Function Definition v1.0
%% Desription
% Uses Horner's algorithm for polynomial evaluation to evaluate the
% curve-fits produced by the CMP-OLS program. The idea behind Horner's rule
% is extended to the evaluation of Laurent polynomials as they may be
% factored identically to regular polynomials. 
%% Future Upgrades
% * Enable Multivariate Polynomial evaluation.
% * Condense numeric evaluation loop.
%% Syntax
% * INPUT(*x*): Discrete data array where to evaluate polynomial.
% * INPUT(*C*): Polynomial coefficients (may have coefficients of negative
% powers). The coefficients must be arranged from lowest order term to
% highest order term.
% * INPUT(*n*): Order of the polynomial. This input is OPTIONAL. If
% unspecified, the function assumes that the order of the polynomial is
% equal to the number of coefficients minus 1. If the specified "n" value
% is less than the number of coefficients provided, the function assumes
% that the extra coefficients are for negative powers.
% * OUTPUT(*y*): Numeric values of the polynomial at the provided x
%% Changelog
% v1.0: Initial Release (Decoupled from "custompolyfit" function).
%% Function Definition
function w = custompolyval(x,C,o,const)
[rX,cX] = size(x);
[rC,cC] = size(C);
[ro,co] = size(o);
if isa(C,'cell') == 0
    error('Coefficient array must be a 1 by Q array of vectors!');
end
if cX ~= cC
    error('Columns of coefficient array do not match those of query array');
end
if cX ~= co
    error('');
end
w = zeros(rX,1);
for i = 1:cC %For all homogeneous polynomials supplied.
    xi = x(1+rX*(i-1):rX*i)';
    mpn = length(C{i}); %Sum of Laurent order and polynomial order.
    m = mpn - o(i); %Laurent order.
    if o(i) > 0 %Compute polynomial component (if any).
        y = ones(rX,1)*C{i}(mpn); % Start with the highest positive power.
        for j = mpn-1:-1:m+1 %Horner's algorithm on the polynomial part.
            y = xi.*y + C{i}(j);
        end
        y = xi.*y;
        w = w + y;
    end
    if m > 0
        y = ones(rX,1)*C{i}(1); %Start with the top-most coefficient.
        for j = 2:m
            y = y./xi + C{i}(j);
        end
        y = y./xi;
        w = w + y;
    end
end
w = w + const;
end