%% Taylor Shift
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE), College of Arts and Sciences (COAS)
%  For use in MA 412, MA 413, AE 435, AE 440 and any other course that
%  would benefit from a data-fitting tool.
%% Description
% Given an array of coefficients of the powers of the indeterminates of a
% polynomial and a desired shift, this routine computes the coefficients of
% the shifted polynomial. This function can be used in applications where
% polynomials of identical shape but differnt position need to be stored
% (like in the construction of NURBS).
%% Formulae
% $f(x) = P_{n}(x)$
%%
% $r(x) = f(Sx)$
%%
% $q(x) = r(x+1)$
%%
% $f(x+S) = q(\frac{x}{S})$
%% Required Plugins
% * binexp.m
%% Changelog
%  v1.0,(08/07/2022): Initial Release. Reinvented the wheel once again.
%% Syntax
% * INPUT(*C_IN*): A numeric row or column array with the coefficients of
% the polynomial sorted by ascending powers.
% * INPUT(*S*): The desired shift. Positive values shift the polynomial in
% the negative direction. Negative values shift the polynomial in the
% positive direction
% * OUTPUT(*C_out*): Coefficients of the shifted polynomial.
%% Function definition
function C_out = tshift(C_in,S)
[r,c] = size(C_in);
n = c+1;%Polynomial's order
if S == 0
    C_out = C_in;
    return;
end
C_out = C_in; %Make a copy of the input coefficients.
pb = 1; %Power buffer.
for ii = 2:c
    pb = pb*S;
    C_out(ii) = pb*C_in(ii);
end

for ii = 2:c
    coeff = binexp(ii-1);
    pb = C_out(ii); %Repurpose the buffer to remember the coefficient.
    for jj = 1:ii-1
        C_out(jj) = C_out(jj) + pb*coeff(jj); 
    end
end
pb = 1; %Reset the power buffer.
for ii = 2:c
    pb = pb/S;
    C_out(ii) = pb*C_out(ii);
end
end