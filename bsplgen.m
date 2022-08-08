%% B-Spline Generator
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE), College of Arts and Sciences (COAS)
%  For use in MA 412, MA 413, AE 435, AE 440 and any other course that
%  would benefit from a data-fitting tool
%% Description
% Given a desired order and  a sufficiently long sequence of nondecreasing
% "knots", computes the piece-wise polynomial coefficients in ascending
% order of a Basis spline (B-spline).
%% Formulae
%  Zeroth-order Basis spline
%%
% $B^{0}_{i}(u) = \begin{array}{c}1\\ 0 \end{array} \quad \begin{array}{c}t_{i} \le u < t_{i+1}\\otherwise\end{array}$
%%
%  Cox-De Boor recursion formula
%%
% $B^{k}_{i}(u) = \left(\frac{u-t_{i}}{t_{i+k}-t_{i}}\right)B^{k-1}_{i}(u) + \left(\frac{t_{i+k+1}-u}{t_{i+k+1}-t_{i+1}}\right)B^{k-1}_{i+1}$
%% Required Plugins
% * horner.m 
%% Changelog
%  v1.0,(08/08/2022): B-spline generation no longer requires the auxiliary
%  buffer.
%  v1.0,(08/07/2022): Initial Release. Reinvented the wheel once again.
%% Syntax
% * INPUT(*kn*): Knot sequence input as either a row or column vector. The
% sequence must be non-decreasing. Repeated knots are allowed as long as
% they and the non-repeated knots are sorted in a non-decreasing order.
% * OUTPUT(*p*): Order of the B-spline. Can be any order. The B-spline
% requires "p+2" knots at minimum. If the user supplies less than "p+2"
% knots, the routine will assume that the first and last knots are repeated
% in a symmetrical way. The repetition depends on the parity of the
% spline's order.
% * OUTPUT(*B*): A numeric, "p+1" by "p+1" array containing the
% coefficients of the B-spline over the range of knots. Each row lists in
% ascending order the coefficients of the piece wise polynomials.
% * OUTPUT(*k_range*): A "p+1" by 2 array that has the knot values that
% define the range of each piece-wise polynomial. The first column lists
% the knot values at the start of the range. The second column lists the
% knot values at the end of the range.
%% Function definition
function [B,k_range,k_idx] = bsplgen(kn,p)
if p < 1 
    error('This routine is for order 1 or above.')
end
[rk,ck] = size(kn);
if rk ~= 1 && ck ~= 1
    error('Must input a row or column array for the knots.');
end
if rk == 1 %Knots input as row vector.
    knots = ck;
end
if ck == 1 %Knots input as column vector.
    knots = rk;
end
if knots < p+2
    error(['Need "p+2" =',num2str(p+2),'knots. Only',num2str(knots),'input']);
end
if knots > p+2
    warning(['Only "p+2" =',num2str(p+2),'knots needed. Detected',num2str(knots-p-2),' in excess. Ignoring.']);
end
for ii = 1:knots-1 %Check if the knot sequence is non-decreasing.
    if kn(ii+1) - kn(ii) < 0
        error('Knot sequence must be non-decreasing.')
    end
end
%Allocate buffer arrays.
B = cell(p+1,1); %Need pointers to buffer arrays.
for ii = 1:p+1
    B{ii} = zeros(p+2-ii,p+2-ii);
    B{ii}(1,1) = 1; %This one is th B^0 spline.
end
%Generate the splines.
for kk = 1:p %Counter "kk" denotes the intermediate order of B-spline being made.
    idx = 1:kk;
    for jj = 1:p-kk+1 %Compute the "k^th" order B-splines into the respective buffers.
        den = kn(jj+kk) - kn(jj); %Denominator
        if den ~= 0 %"Bootstrap" component.
            B{jj}(idx,idx) =  -B{jj}(idx,idx)*kn(jj)/den;
            B{jj}(idx,idx+1) = B{jj}(idx,idx+1) - B{jj}(idx,idx)/kn(jj);
        end
        den = kn(jj+kk+1) - kn(jj+1); %Denominator
        if den ~= 0 %"Next" spline component.
            B{jj}(idx+1,idx) = B{jj}(idx+1,idx) + B{jj+1}(idx,idx)*kn(jj+1+kk)/den;
            B{jj}(idx+1,idx+1) = B{jj}(idx+1,idx+1) - B{jj+1}(idx,idx)/den;
        end
    end
end
if nargout > 1
    k_range = zeros(knots-1,2);
    for kk = 1:knots-1
        k_range(kk,1) = kn(kk);
        k_range(kk,2) = kn(kk+1);
    end
end
if nargout > 2
    k_idx = zeros(knots-1,2);
    for kk = 1:knots-1
        k_idx(kk,1) = kk;
        k_idx(kk,2) = kk + 1;
    end
end
end