%% Cubic Spline Coefficient Generator v1.1
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%% Description
% A function that takes in independent (x) and dependent data (y) sets with
% $n$ observations and computes the associated $4(n-1)$ polynomial
% coefficients of the associated $n-1$ splines. "splinegen.m" outputs the
% coefficients as a $n-1 \times 4$ array for each combination of x and y
% data sets. An optional parameter gives the user control over the
% 2nd-derivative values at the data's endpoints to create non-natural
% splines if necessary. Another optional parameter further outputs the
% coefficients of the corresponding 1st-derivatives of the splines.
%% Required Plugins
% * tridiag.m (A tridiagonal system needs to be solved to generate splines)
%% Changelog
%  v1.1,(06/13/2021): Added saving functionality. Can export spline
%  coefficients as several .csv files or as a multi-sheet .xlsx file.
%  Customization of naming convention in the form of y = S(x) "y equals
%  spline of x," where "y" and "x" are the variable names which the user
%  may customize.
%  v1.0,(06/01/2021): Initial Release.
%% Syntax
% * INPUT(*X*): Independent data set(s), each of size $N \times 1$,
% concatenated horizontally.
% * INPUT(*Y*): Dependent data set(s), each of size $N \times 1$
% concatenated horizontally.
% * INPUT(*varargin*): Variable arguments of the form "key, value:"
% * KEY/VALUE ('filename', *filename*): The string 'filename' followed by a
% user-entered string that denotes the save filename. *filename* can
% contain file extensions to specify the saving format (i.e. .csv or .xlsx)
% If the *filename* does not include an extension, the default is .csv.
% * KEV/VALUE ('Xname', *xname*): A cell array that contains custom
% independent data names. User must be careful to ensure that there are as
% many custom names as there are indepedendent data sets. *NOTE:* this
% input is redundant if *filename* is unspecified.
% * KEY/VALUE ('Yname', *yname*): Same as *xname* but for the dependent data
% sets.
% * KEY/VALUE ('clamps', *clamps*): A two-element numeric vector that
% specifies the values of the second derivative at the endpoints so that
% the user may create "non-natural" splines.
% * OUTPUT(*fc*): Cell array whose elements contain the spline
% coefficients as $N-1 \times 4$ matrices.
%% Function Definition
function fc = splinegen(X,Y,varargin)
%Dummy variables.
clamps = [0,0]; %Assume Natural splines by default.
xname = 0; %Dummy variable for X-names.
yname = 0; %Dummy variable for Y-names.
filename = 0; %Dummy variable for filename.
format = '.csv'; %Default file extension is comma-separated-value(s).

[Nx,nx] = size(X); %Observations and # of data sets of independent-variable.
[Ny,ny] = size(Y); %Observations and # of data sets of dependent-variable.
if Ny ~= Nx
    error('Unequal observations! (rows of x do not match those of y)')
end

%Variable Argument recognition and error checking.
inpairs = mod(nargin,2); %"Incomplete" pair.
if inpairs == 1
    error('Uneven number of inputs!')
end
pairs = nargin/2; %Number of "key-value" argument pairs.
if pairs > 6
    error('Excessive number inputs! This function deals with only five(5).')
elseif pairs < 1
    error('Insufficient inputs. At minimum must input X and Y data!')
end
for i = 1:pairs-1 %For each argument pair.
    key = 2*i-1; %index of character array denoting a "key"
    val = 2*i; %index of character array denoting key's value.
    switch varargin{key} %Attempt to recognize valid inputs.
        case 'clamps' %If input recognized as custom endpoint "S" values.
            clamps = varargin{val};
            if isnumeric(clamps) ~= 1
                error('Input for "clamps" must be a numeric 1 by 2 array.')
            end
        case {'Xname','Yname'} %If input recognized as names of independent data set(s).
            if varargin{key} == 'Xname'
                xname = varargin{val};
                lim = nx;
                str = 'X';
            else
                yname = varargin{val};
                lim = ny;
                str = 'Y';
            end
            if iscell(varargin{val}) == 1 %If the input is a cell.
                if length(varargin{val}) == lim %Number of cell entries must match that of corresponding data sets.
                    for j  = 1:length(varargin{val}) %For each cell entry...
                        if ischar(varargin{val}{j}) == 0 %check that ALL are character arrays
                            error('At least one name entry is not a character array!')
                        end
                    end
                else %Throw mismatch error
                    error(['rows of ',varargin{key},' must match those of ',str])
                end
            else
                error([varargin{key},'input must be of class "cell" with all "char" entries'])
            end
        case 'filename' %If input recognized as filename of save.
            filename = varargin{val};
            if ischar(filename) ~= 1 %Filename must be input as a character array.
                error('filename must be a character array!')
            else
                L = length(filename); %Get filename's number of characters.
                k = 0; %A counter used to get the index of the file extension.
                while k < L
                    if filename(L-k) == '.'
                        if sum(strcmpi(filename(L-k:L),{'.csv','.xlsx'})) == 0
                            warning('Unrecognizable format extension. Output will save as ".csv"')
                        else
                            format = filename(L-k:L); %Overwrite the default .csv with the detected format.
                        end
                        filename(L-k:L) = []; %Erase the extension from the filename.
                        k = L + 1; %Exit the while loop.
                    else
                        k = k + 1; %Update filename letter counter.
                    end
                end
                if k == L
                    warning('No file extension detected in filename. Output(s) will save as ".csv"')
                end
            end
        otherwise
            error('Unrecognizable input. Valid options: "filename", "Xname", and "Yname".')
    end
end

%Memory allocation.
fc = cell(ny,nx); %Preallocate memory for output.
b = zeros(Nx,1); %Preallocate memory for the RHS vector.
A = zeros(Nx,3); %Preallocate memory for A in sparse mode.
if filename ~= 0 %If a filename was specified, define default output names.
    if iscell(xname) == 0 %No X-variable names specified
        xname = cell(nx,1);
        for i = 1:nx
            xname{i} = ['X',num2str(i)];
        end
    end
    if iscell(yname) == 0 %No Y-variable names specified
        yname = cell(ny,1);
        for j = 1:ny
            yname{j} = ['Y',num2str(j)];
        end
    end
end

%Additional pre-loop computations.
A([Nx+1,2*Nx]) = [1,1]; %Coefficient of the endpoints (sparse mode)
b([1,Nx]) = clamps; %Value of the left end clamp (Natural splines). [CURRENTLY REDUNDANT]
internal = (2:Nx-1)'; %LI of internal RHS elements (endpoints excluded).
backward = internal - 1; %LI of RHS elements excluding last two.
forward = internal + 1; %LI of RHS elements excluding first two.
mainint = internal + Nx; %LI of internal A main diagonal elements.
supint = mainint + Nx; %LI of internal A superdiagonal elements.
backinc = (1:Nx-1)'; %LI of RHS elements excluding last one.
forwinc = backinc + 1; %LI of RHS elements excluding first one.

for i = 1:nx %For all independent data sets.
    %Independent data vectorizations.
    xpacer = Nx*(i-1); %LI to discriminate independent data sets.
    xim1 = backward + xpacer; %"i-1" data point vectorization (excludes last two).
    xi = internal + xpacer; %"i" data point vectorization (excludes endpoints).
    xip1 = forward + xpacer; %"i+1" data point vectorization (excludes first two).
    xiinc = backinc + xpacer; %"i" data point vectorization (excludes only right end point).
    xip1inc = xiinc + 1; %"i+1" data point vectorization (excludes only left end point).
    
    A(internal) = X(xi) - X(xim1); %Update internal subdiagonal coefficients.
    A(mainint) = 2*(X(xip1) - X(xim1)); %Update internal main diagonal coefficients.
    A(supint) = X(xip1) - X(xi); %Update internal superdiagonal coefficients.
    deltaX = X(xip1inc) - X(xiinc);
    for j = 1:ny %Generate spline coefficients for all dependent data sets.
        %Independent data vectorizations.
        ypacer = Ny*(j-1); %LI to discriminate dependent data sets.
        yim1 = backward + ypacer; %Same as above but for y.
        yi = internal + ypacer; %Same as above but for y.
        yip1 = forward + ypacer; %Same as above but for y.
        yiinc = backinc + ypacer; %Same as above but for y.
        yip1inc = yiinc + 1; %Same as above but for y.
        
        %RHS editing and tridiagonal system solution.
        b(internal) = 6*((Y(yip1) - Y(yi))./A(supint)...
            - (Y(yi) - Y(yim1))./A(internal)); %Update RHS vector.
        S = tridiag(A,b); %Callout to Thomas algorithm solver.
        
        %Vectorized Cubic Spline coefficients (Yp = k3*x^3 + k2*x^2 + k1*x + k0).
        k3 = (S(forwinc) - S(backinc))./(6*deltaX); %3rd-order coefficients.
        k2 = (X(xip1inc).*S(backinc) - X(xiinc).*S(forwinc))./(2*deltaX); %2nd-order coefficients.
        k1 = (Y(yip1inc) - Y(yiinc))./deltaX...
            - k3.*(X(xip1inc).^2 + X(xip1inc).*X(xiinc) + X(xiinc).^2)...
            - k2.*(X(xip1inc)+ X(xiinc)); %1st-order coefficients.
        k0 = Y(yiinc) - X(xiinc).*(k1 + X(xiinc).*(k2 + X(xiinc).*k3)); %constant coefficients
        dataset = j + ny*(i-1); %LI of current set of splines within cell.
        fc{dataset} = [k3,k2,k1,k0]; %Copy results to output cell array.

        %Output data on the go.
        if filename ~= 0
            switch format
                case '.csv' %If user specified comma-separated-values
                    csvwrite([filename,'_',yname{j},'= S(',xname{i},')',format],fc{dataset})
                case '.xlsx'
                    xlswrite(filename,fc{dataset},[yname{j},'= S(',xname{i},')'])
            end
        end
    end
end
end