%% Multi-Dimensional, Multi-guess Newton-Raphson (MDMGNR)Function Definition v2.0
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%% Description
% A function that implements the exact Newton method for multiple
% dimensions and guesses. Requires the system of nonlinear equations and
% its Jacobian matrix to be input as function handles. A matrix inversion
% algorithm is required for solving multidimensional problems.
%% Required Plugins
% * Gauss Jordan matrix inversion algorithm.
%% Changelog
%  v1.2,(11/15/2020): Added optional "par" input to reduce overhead from
%  function handles that define constant variables in their respective
%  scripts.
%  v1.1,(11/11/2020): Enabled multidimensional root finding.
%  v1.0,(10/02/2020): Initial Release (unidimensional root-finding only).
%% Syntax
% * INPUT(*fval*): Vector of values that fvec (below) equals to.
% * INPUT(*fvec*): Function handle vector (or scalar). NOTE: f(x,y,z) = fvec(x,y,z) - fval
% * INPUT(*fmat*): Jacobian matrix (or scalar derivative) function handle.
% * INPUT(*Xn*): Initial guess vector (can supply multiple initial guesses
% as a row vector.
% * INPUT(*par*): Optional aditional parameters used by fmat and fvec (must
% be a structured array).
% * INPUT(*tol*): Maximum error (based on successive iterations). Computed
% as the the difference of successive iterations.
% * INPUT(*maxit*): Maximum allowed iterations
% * OUTPUT(*Xnp1*): Root prediction vector.
% * OUTPUT(*k*): Vector of iterations at convergence per prediction.
% * OUTPUT(*delta*): Error accrued in prediction vector at convergence.
%% Formulae
% * $$X_{n+1} = X_{n} - \frac{f(X_{n})}{f^{\prime}(X_{n})}$
% * $$\vec{X}_{n+1} = \vec{X}_{n} - J^{-1} \times \vec{f}(\vec{X}_{n})$
%% Function Definition

%something = newt('F',@AoAs,'X0',[1;2;3])
function [Xnp1,k,res] = newt(varargin)
%Input registers
tol_spec = 0;%tolerance register.
maxit_spec = 0;%maximum iterations register.
par_spec = 0;%optional parameters register.
F_spec = 0;%Function vector (F(xk)) register.
J_spec = 0;%Jacobian matrix (J(xk)) register.
X0_spec = 0;%Initial guess vector (x0) register.
conv_spec = 0;%Convergence criteria register.
norm_spec = 0;%Vector norm register.
hist_spec = 0;%Convergence history register.

%Input collection and error checking.
if mod(nargin,2) == 1 %Incomplete input pair.
    error('Uneven number of inputs! Must supply complete pairs.')
elseif nargin > 12 %Excessive number of inputs.
    error('Excessive number of input pairs! This function deals with only five(5).')
elseif nargin < 2 %Insufficient number of inputs.
    error('Insufficient input pairs! (minimum of two(2) physical inputs!)')
else %Viable number of input pairs given.
    pairs = nargin/2; %Number of "key-value" argument pairs.
    for i = 1:pairs %For each argument pair.
        switch varargin{2*i-1} %Attempt to recognize valid inputs.
            case 'tol'%tolerance on prescribed norm.
                if isnumeric(varargin{2*i}) == 1
                    [r,c] = size(varargin{2*i});%Check the size of the input.
                    if r == c && c == 1%Ensure the numeric input is a 1by1 array.
                        tol = varargin{2*i};
                        tol_spec = 1;%Register that tol was input.
                    else
                        error('ERROR: "tol" must be scalar, not an array!')
                    end
                else
                    error('ERROR: "tol" input must be numeric!')
                end
            case 'norm'%norm type to establish convergence
                if isnumeric(varargin{2*i})==1
                    norm_spec = 1;
                    switch varargin{2*i}
                        case 1
                            norm = @(V) (sum(abs(V)));%Taxicab norm (p=1)
                        case 2
                            norm = @(V) (sqrt(sum(V.*V)));%Euclidean norm (p=2)
                        case Inf
                            norm = @(V) (max(abs(V)));%Infinity norm (p=Inf)
                        otherwise
                            error('ERROR: Invalid "norm" input. Valid options are "1", "2", and "Inf".')
                    end
                else
                    error('ERROR: "norm" must be input as numeric 1, 2,, or Inf!')
                end
            case 'conv' %Convergence criteria (i.e., how compare norms againts tolerance)
                switch varargin{2*i}
                    case 'a'%absolute norm on the iterate.
                    case 'r'%relative norm on the iterate.
                    case 's'%relative norm on the step.
                    otherwise
                        error('ERROR: Invalid "norm" input. Valid options are "r", "a", and "s".')
                end
                    conv_spec=1;
            case 'maxit'%maximum newton iterations
                if isnumeric(varargin{2*i}) == 1%Check to see if the input is numeric.
                    [r,c] = size(varargin{2*i});%Check the size of the input.
                    if r == c && c == 1%Ensure the numeric input is a 1by1 array.
                        if mod(varargin{2*i},1)==0%Ensure the scalar is an integer.
                            maxit = varargin{2*i};%Collect user input.
                            maxit_spec = 1;%Register that maxit was input.
                        else
                            error('ERROR: Input scalar is not an integer!')
                        end
                    else
                        error('ERROR: "maxit" must be an integer scalar, not an array!')
                    end
                else
                    error('ERROR: "maxit" input must be numeric (specifically an integer)!')
                end
            case 'par'%Optional data
                if isstruct(varargin{2*i}) == 1
                    par_spec = 1;
                    par = varargin{2*i};
                else
                    error('ERROR: Optional data must be passed as a structured array!');
                end
            case 'F'%Function
                if ishandle(varargin{2*i})
                    F_spec = 1;
                    F = varargin{2*i};
                else
                    error('ERROR: Function vector must be specified as a handle!')
                end
            case 'J'%Jacobian
                if ishandle(varargin{2*i})
                    J_spec = 1;
                    J = varargin{2*i};
                else
                    error('ERROR: Jacobian must be specified as a handle!')
                end
            case 'x0'
                X0_spec = 1;
                xn = varargin{2*i};
            case 'fval'
            case 'hist'%Convergence history
                hist_spec = 1;
        end
    end
end
%Default values.
if F_spec == 0
    error('ERROR: No function vector F(x_k) specified! (Nothing to solve');
end
if X0_spec == 0
    Xn = 0;
end
if J_spec == 0
    %J = @(F,Xk) (for end );%Set up a finite difference scheme.
end
%Program defaults
if tol_spec == 0%If tolerance was unspecified.
    tol = 1e-6;
end
if maxit_spec == 0%If maximum number of iterations was unspecified.
    maxit = 100;
end
if par_spec == 0
    par = struct('y',1.4);
end
if norm_spec == 0%If norm type was unspecified.
    norm = @(V) (sqrt(sum(V.*V)));%Euclidean norm (p=2)
end









[rXn,cXn] = size(Xn);
%cXn denotes the number of initial guesses.
%rXn denotes the number of unknowns to solve.
Xnp1 = zeros(rXn,cXn); %Preallocate memory for the Output matrix.
k = zeros(1,cXn); %Preallocate memory for the iteration counters.
res = ones(1,cXn); %Preallocate memory for the error trackers.
for d = 1:cXn
    while res(d) > tol && k(d) <= maxit %Enforce convergence criteria.
        Fx = F(Xn(:,d));
        s = J(Xn(:,d))\(-Fx);%Compute the newton step.
        Xnp1(:,d) = s + Xn(:,d);
        res = conv;
        Xn(:,d) = Xnp1(:,d);
        k(d) = k(d)+1;%Update iteration counter.
    end
end
    function jacobian = fd(F,Fx,X)%finite difference function definition.
        [a,~] = size(X);%Get dimensions of input vector.
        delta = sqrt(eps);%FD step is square root of machine epsilon.
        jacobian = zeros(a,a);
        for j = 1:a
            X(j) = X(j)+delta;
            jacobian((j:a)+a*(j-1)) = (F(X)-Fx)/delta;
            X(j) = X(j)-delta;
        end
    end
end