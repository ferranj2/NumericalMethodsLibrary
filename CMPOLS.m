%% Constrained, Multidimesional, Polynomial Ordinary Least Squares Regression (CMP-OLS) v1.0
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE), College of Arts and Sciences (COAS)
%  For use in MA 412, MA 413, AE 435, AE 440 and any other course that
%  would benefit from a data-fitting tool.
%% Description
% This tool (an algorithm really) was made to assist students and faculty
% who exhibit the need to approximate as a continuous function, a quantity
% or property that is originally only available as discrete observations.
% The purpose of this tool is two-fold: Firstly, to provide users with a
% simple and versatile mechanism to curve-fit discrete data, and secondly,
% to approximate the inverse of transcendental or nonlinear functions so as
% to generate accurate initial guess schema for use in nonlinear
% Newton-Rhapson solvers. The program is versatile because it is
% multidimensional in nature, offers the user to impose equality
% constraints on their curve fit via Lagrange multipliers, and because it
% can incorporate negative powers as regressors. The brunt of this
% curve-fitting tool is a custom-made algorithm for building the design
% matrix $X^{T}X$ associated with OLS regression. The algorithm seeks to
% produce $X^{T}X$ by leveraging as many of its features that result from
% it being a symmetric tensor.
%% Syntax
% * INPUT(*x*): Independent variable observation array. This contains the
% discrete data of the explanatory variables. The array's rows denote the
% observations whereas its columns denote additional explanatory variables.
% * INPUT(*y*): Dependent variable observation array. Contains the discrete
% data for the reponse variable.
% * INPUT(*o*): Polynomial/Laurent order matrix. This specifies, for each
% predictor variable, the highest positive and negative powers to include
% in the OLS regression.
% * INPUT(*c*): Constraint array. This array contains coordinates of the
% n-dimensional points that MUST be included in the curve-fit. The
% resulting curve-fit will perfectly (as far as computational arithmetic
% precision will permit) go through these points.
% * OUTPUT(*coeff*): A vector array that contains the coefficients of the
% regressors. These should be saved so that the curve-fit may be evaluated
% at later times.
%% Future Upgrades
% * Enable equality constraints on the gradient of the regression field.
% * Enable regression with coupled terms (e.g.  bilinears x*y, x*z, y*z;
% biquadratics x^2*y, x^2*z, x^2*y^2, x^2*z^2, etc.
% * Incorporate statistical analysis from ANOVA and ANCOVA.
% * Include data visualization (MATLAB plotting).
clc
clear
close all
%% Test Inputs
x1 = [1,5,3,7,6,2,6,8,9,1,1,2,4]'; %sum = 55
x2 = [3,5,7,9,1,6,6,8,1,4,2,9,7]'; %sum = 68
x3 = [1,9,3,9,5,9,2,9,8,3,5,8,4]'; %sum = 75
Y = [1;2;3;4;5;6;7;8;9;10;11;12;13]; %sum = 91
X = [x1,x2,x3];
o = [3,0,3;...
    2,2,0];
K = [1,2,3;...
    2,3,4];
K = [];
%function CMPOLS
[Nx,Qx] = size(X);%Number of predictor observations (Nx) & variables (Q).
[Ny,Qy] = size(Y);%Number of dependent observations (Ny)
[Nk,Qk] = size(K);%Number of constraints imposed on the system.
[No,Qo] = size(o);
sPn = sum(o(1:2:2*Qx-1));%Sum of polynomial orders of all x_i.
sLn = sum(o(2:2:2*Qx));%Sum of Laurent orders of all x_i.
dof = Ny - sLn - sPn - 1 - Nk;%Degrees of freedom.
%% Multivariate XTX builder (bootstrap) (Generate massive rectangular arrays)
%{
R = sLn + 1 + sPn; %Number of unkown coefficients.
P = 1 + sLn*(1 + R); %Single-Index location of the Pivotal coefficient.
BnP = P - sLn*(R - 1); %Pivotal C-band number.
A = zeros(R,R); %"Blank Canvas" Vandermonde matrix.
RHS = zeros(R,1); %"Blank Canvas" the RHS vector.
generatrix = ones(Nx,sLn + sPn); %"Blank Canvas" an array of sums.
Pk = ((1:Nx) + Nx*(sLn-1))'; %Poly-counter.
Lk = ((1:Nx) + Nx*(sLn))'; %Laur-counter.
vec_temp = zeros(Nx,1); %"Blank Canvas" a "Dummy" vector of sums.
%Generatrix builder
for i = 1:Q %For all independent variables.
    p_i = 1 + 2*(i-1); %Linear index of current variable's polynomial order.
    x_i = ((1:Nx) + Nx*(i-1))'; %Linear index of the observations of the "ith" independent variable.
    Pni = n(p_i); %Extract current Variable's polynomial order.
    Lni = n(p_i + 1); %Extract current variable's Laurent order.
    if Pni > Lni
        hmp = Pni; %Polynomial order is the "higher magnitude power."
    else
        hmp = Lni; %Laurent order is the "higher magnitude power."
    end
    vec_temp(:) = 1; %Reset dummy vector for next variable.
    for j = 1:hmp
        vec_temp = vec_temp.*x(x_i); %Increase the order of the dummy vector.
        if Pni > 0 && j <= Pni
            Pk = Pk + Nx; %Update poly-counter.
            generatrix(Pk) = vec_temp; %Populate polynomial columns.
        end
        if Lni > 0 && j <= Lni
            Lk = Lk - Nx; %Update laur-counter.
            generatrix(Lk) = 1./vec_temp; %Populate laurent columns.
        end
    end
end
generatrix
%}
%% Error Checking
if Nx ~= Ny
    error('Number of independent observations does not equal that of dependent observations (Nx ~= Ny).');
end
if Qo ~= Qx
    error('Dimensionality of array of orders does not match that of predictors.');
end
if No >2
    error('Only polynomial and Laurent orders supported.')
end
if Nk > 0 && Qk ~= Qx
    error('Dimensionality of constraints imposed does not match that of predictors.');
end
if dof == 0
    warning('dof = 0! Regression analysis impossible, direct fit will ensue.');
elseif dof < 0
    error('Underdetermined system (dof < 0).')
end
%Preallocate memory for the output structured arrays.
describe = struct(...
    'min',zeros(1,1+Qx),...%Minimum values
    'Q1',zeros(1,1+Qx),... %First quartile
    'med',zeros(1,1+Qx),...%Median (second quartile)
    'Q3',zeros(1,1+Qx),... %Third quartile
    'max',zeros(1,1+Qx),...%Maximum values
    'avg',zeros(1,1+Qx),...%Arithmetic means.
    'std',zeros(1,1+Qx),...%Standard deviation.
    'var',zeros(1,1+Qx),...%Variance.
    'skw',zeros(1,1+Qx),...%Skewness.
    'kur',zeros(1,1+Qx));  %Kurtosis.
data = struct(...
    'X',X,...              %Input predictor data.
    'Y',Y,...              %Input dependent data.
    'Zx',zeros(Ny,Qx),...  %Standardized predictor data.
    'Zy',zeros(Ny,Qy),...  %Standardized dependent data.
    'res',zeros(Ny,Qy));     %Residuals.
OLS = struct(...
    'SSR',0,...            %Sum of Squares of residuals.
    'SSE',0,...            %Sum of Squares of errors.
    'R2',0,...             %Coefficient of Determination.
    'R2A',0,...            %Adjusted R2 value.
    'N',0,...              %Number of observations.
    'K',0,...              %Number of constraints.
    'dof',0,...            %Degrees of freedom.
    'con',0,...            %Condition number of the XTX matrix.
    'C',{cell(1,Qx)});       %Coefficients of the regressors.
for i = 1:Qx+1
    if i == 1
        array = Y;
    else
        array = X(1+Ny*(i-2):Ny*(i-1));
    end
    describe.min(i) = min(array);
    %describe.Q1(i) = ;
    describe.med(i) = median(array);
    %describe.Q3(i) = ;
    describe.max(i) = max(array);
    describe.avg(i) = sum(array)/Ny;
    describe.var(i) = sum((array-describe.avg(i)).^2)/Ny;
    describe.std(i) = sqrt(describe.var(i));
    %describe.skw(i) = sum(((x(1+Ny(i-1):Ny*i)-describe.avg(i))./describe.std(i)).^3)/Ny;
    describe.skw(i) = sum((array-describe.avg(i)).^3)/(Ny*describe.var(i)^(3/2));
    describe.kur(i) = sum((array-describe.avg(i)).^4)/(Ny*describe.var(i)^2);
end
%% Multivariate XTX builder (originator) (Pivotal extrusion + transposition)
R = sLn + 1 + sPn + Nk; %Size of XTX.
Rp1 = R + 1;%Takes a diagonal step.
Rm1 = R - 1;%Takes a counter diagonal step.
XTX = zeros(R,R); %Preallocate memory for the design matrix.
RHS = zeros(R,1); %Preallocate memory for the RHS vector.
N = 0; %Counter for polynomial orders spanned.
M = 0; %Counter for Laurent orders spanned.
if Qx > 1 && sPn > 0
    P_el = zeros(Ny,sPn); %Preallocate memory for positive primordial lists.
end
if Qx > 1 && sLn > 0
    L_el = zeros(Ny,sLn); %Preallocate memory for negative primordial lists.
end
temp = zeros(Nx,1);%Preallocate memory for a vector of powers.
P = 1 + sLn*Rp1;%Linear index of the Pivotal coefficient.
RHS(sLn+1) = sum(Y);
XTX(P) = Nx;%Pivot element (equals Nx always).
for i = 1:Qx %For each independent variable.
    p_i = 1 + 2*(i-1);%LI of current variable's polynomial order.
    n = o(p_i);%Current variable's polynomial order.
    m = o(p_i+1);%Current variable's laurent order.
    N = N + n;%Update polynomial orders spanned.
    M = M + m;%Update Laurent orders spanned.
    %TL = Top left
    %TR = Top right
    %BL = Bottom left
    %BR = Bottom right
    %sp = start point
    %ep = end point
    %hp = high power
    %lp = low power.
    %cbandn = counter-diagonal band with positive power sums.
    %cbandm = counter-diagonal band with negative power sums.
    temp(:) = 1; %Reset vector of powers.
    xi = (1+Ny*(i-1):Ny*i)'; %LI of observations of ith independent variable.
    
    if m == 0 && n == 0
        error('An explanatory variable has m = n = 0!');
    elseif m == 0 && n > 0 %If variable has a polynomial order only.
        TL = P; %Top-Left corner coincides with the pivotal coefficient.
        BL = TL + N;%BL below TL by polynomial order spanned.
        TR = TL + R*N;%
        BR = BL + R*N;
        sp = [TL,... %corner.
            BL-n+1:BL,... %Left vertical edge (+jump).
            BR-R*(n-1):R:BR]; %Bottom horizontal edge (no jump).
        ep = [TL,... %corner.
            TR-R*(n-1):R:TR,... %Top horizontal edge.
            BR-n+1:BR]; %Right vertical edge.
        for j = 2:2*n+1 %Skip pivot as it has already been filled in.
            temp = temp.*X(xi);
            if j == 2
                cbandn = [sp(j),ep(j)];
            elseif j <= n+1
                cbandn = [sp(j),...
                    sp(j)+(N-n+1)*R-1:Rm1:ep(j)+N-n+1-R,...
                    ep(j)];
            elseif j > n+1
                cbandn = sp(j):Rm1:ep(j);
            end
            XTX(cbandn) = sum(temp);
            if j <= n+1 && Qx > 1 && sPn > 0
                P_el(1+Ny*(N-n+j-2):Ny*(N-n+j-1)) = temp;%Store basis list.
                RHS(sLn+1+N-n+j-1) = sum(Y.*temp);
            end
        end
        if i > 1%Compute mixed term sums from the basis lists.
            if N > 0 %If polynomial orders spanned.
                for j = N-n+1:N
                    for k = 1:N-n
                        mix = P+j+R*k;
                        XTX(mix) = sum(P_el(1+Ny*(j-1):Ny*j).*P_el(1+Ny*(k-1):Ny*k));
                        XTX(P+R*j+k) = XTX(mix);
                    end
                end
            end
            if M > 0 %If Laurent orders spanned.
                for j = N-n+1:N
                    for k = 1:M
                        mix = P+j-R*k;
                        XTX(mix) = sum(P_el(1+Ny*(j-1):Ny*j).*L_el(1+Ny*(k-1):Ny*k));
                        XTX(P+R*j-k) = XTX(mix);
                    end
                end
            end
        end
    elseif n == 0 && m > 0 %If variable has a Laurent order only.
        TL = P - Rp1*M;
        BL = TL + M;
        TR = TL + R*M;
        BR = P; %Bottom-Right corner is the pivotal coefficient.
        sp = [TL:TL+m-1,... %Left vertical edge.
            BL:R:BL+R*(m-1),... %Bottom horizontal edge.
            BR]; %corner
        ep = [TL:R:TL+R*(m-1),... %Top horizontal edge.
            TR:TR+m-1,... %Right vertical edge.
            BR]; %corner
        for j = 2*m:-1:1 %Skip pivot as it has already been filled in.
            temp = temp./X(xi); %Update temporary vector.
            if j == 2*m
                cbandm = [sp(j),ep(j)];
            elseif j >= m
                cbandm = [sp(j),...
                    sp(j)-(M-m+1)+R:Rm1:ep(j)-R*(M-m+1)+1,...
                    ep(j)];
            elseif j < m
                cbandm =sp(j):Rm1:ep(j);
            end
            XTX(cbandm) = sum(temp);
            if j > m && Qx > 1 && sLn > 0
                L_el(1+Ny*(M-m-j+2*m):Ny*(M-m-j+2*m+1)) = temp;%Store basis list.
                RHS(sLn+1-M+1+j-2*m) = sum(Y.*temp);
            end
        end
        if i > 1%Compute mixed term sums from the primordial lists.
            if N > 0%If polynomial orders were spanned.
                for j = M-m+1:M
                    for k = 1:N
                        mix = P+k-R*j;
                        XTX(mix) = sum(L_el(1+Ny*(j-1):Ny*j).*P_el(1+Ny*(k-1):Ny*k));
                        XTX(P+R*k-j) = XTX(mix);
                    end
                end
            end
            if M > 0%If Laurent orders were spanned.
                for j = M-m+1:M
                    for k = 1:M-m
                        mix = P-k-R*j;
                        XTX(mix) = sum(L_el(1+Ny*(j-1):Ny*j).*L_el(1+Ny*(k-1):Ny*k));
                        XTX(P-k*R-j) = XTX(mix);
                    end
                end
            end
        end
    elseif n > 0 && m > 0 %If variable has both polynomial and Laurent order.
        TL = P-Rp1*M;
        BL = TL+M+N;
        TR = TL+R*(M+N);
        BR = P+Rp1*N;
        sp = [TL,... %Top-Left Corner.
            TL+1:TL+m-1,... %Top-Left vertical segment.
            TL+M,... %Left edge intersects the pivotal row.
            BL-n+1:BL-1,... %Bottom-Left vertical segment.
            BL,... %Bottom-Left Corner.
            BL+R:R:BL+R*(m-1),... %Bottom-Left horizontal segment.
            BL+R*M,... %Bottom edge intersects the pivotal column.
            BR-R*(n-1):R:BR-R,... %Bottom-Right horizontal segment.
            BR]; %Bottom-Right corner.
        ep = [TL,...%Top-Left Corner.
            TL+R:R:TL+R*(m-1),... %Top-Left horizontal segment.
            TL+R*M,... %Top edge intersects pivotal column.
            TR-R*(n-1):R:TR-R,...%Top-Right segment.
            TR,... %Top-Right corner.
            TR+1:TR+m-1,... %Top-Right vertical segment.
            TR + M,... %Right edge intersects pivotal row.
            BR-n+1:BR-1,... %Bottom-Right vertical segment.
            BR]; %Bottom-Right corner.
        if n > m %if polynomial order is greater than Laurent one.
            hp = n; %n is the high power.
            lp = m; %m is the low power.
        elseif m > n %If Laurent order greater than polynomial one.
            hp = m; %m is the high power.
            lp = n; %n is the low power.
        else %m and n are equal, both are the high and low powers.
            hp = n;
            lp = n;
        end
        BnP = 1+2*m; %"Local" pivotal band number.
        for j = 1:2*hp
            temp = temp.*X(xi);
            if j<=2*n
                Bn = BnP+j;%Band number of positive power sums.
                if j == 1
                    cbandn = [sp(Bn):Rm1:sp(Bn)+Rm1*(n-j-1),...
                        P+N-n+j,...
                        P+(N-n+j)*R,...
                        ep(Bn)-Rm1*(m-j):Rm1:ep(Bn)];
                elseif n-j >= 0
                    cbandn =[sp(Bn):Rm1:sp(Bn)+Rm1*(n-j-1),...
                        P+N-n+j,...
                        P+N-n+j+(N-n+1)*R-1:Rm1:P+(N-n+j)*R+(N-n+1)-R,...
                        P+(N-n+j)*R,...
                        P+(N-n+j)*R-(M-m+1)+R:Rm1:ep(Bn)];
                elseif j == n
                    cbandn = [sp(Bn),...
                        sp(Bn)+(N-n+1)*R-1:Rm1:ep(Bn)+(N-n+1)-R,...
                        ep(Bn)];
                else
                    cbandn = sp(Bn):Rm1:ep(Bn);
                end
            end
            XTX(cbandn) = sum(temp);%Populate c-band of posite power sums.
            if j<=n && Qx > 1
                P_el(1+Ny*(N-n+j-1):Ny*(N-n+j)) = temp;%Store primordial list.
                RHS(sLn+1+N-n+j) = sum(Y.*temp);
            end
            if j <= 2*m
                Bn = BnP-j;%Band numbers for negative power sums.
                if j == 1
                    cbandm = [sp(Bn):Rm1:sp(Bn)+Rm1*(m-j-1),...
                        P-(M-m+j)*R,...
                        P-(M-m+j),...
                        ep(Bn)-Rm1*(m-j-1):Rm1:ep(Bn)];
                elseif j <= m
                    cbandm = [sp(Bn):Rm1:sp(Bn)+Rm1*(m-j-1),...
                        P-(M-m+j)*R,...
                        P-(M-m+j)*R-(M-m+1)+R:Rm1:P-(M-m+j)-(M-m+1)*R+1,...
                        P-(M-m+j),...
                        P-(M-m+j)+R*(N-n+1)-1:Rm1:ep(Bn)];
                elseif j > m
                    cbandm = sp(Bn):Rm1:ep(Bn);
                end
            end
            %{
            if i ==3
                sum(1./(x(xi).^j))
                sum(1./temp)
                %[1./x(xi).^j,1./temp]
            end
            %}
            XTX(cbandm) = sum(1./temp);
            if j <= m && Qx > 1
                L_el(1+Ny*(M-m+j-1):Ny*(M-m+j)) = 1./temp;%Store base list.
                RHS(sLn+1-M+m-j) = sum(Y./temp);
            end
        end
        XTX([sp(BnP):Rm1:sp(BnP)+Rm1*(lp-1),ep(BnP)-Rm1*(lp-1):Rm1:ep(BnP)]) = Ny;%Populate local pivotal band.
        if i > 1
            if N > 0
                for j = N-n+1:N
                    %Mixed sums of Polynomial base with other polynomial
                    %bases.
                    for k = 1:N-n
                        mix = P + j + R*k;
                        XTX(mix) = sum(P_el(1+Ny*(j-1):Ny*j).*P_el(1+Ny*(k-1):Ny*k));
                        XTX(P+R*j+k) = XTX(mix);
                    end
                    %Mixed sums of Polynomial base with other Laurent
                    %bases.
                    for k = 1:M-m
                        mix = P + j - R*k;
                        XTX(mix) = sum(P_el(1+Ny*(j-1):Ny*j).*L_el(1+Ny*(k-1):Ny*k));
                        XTX(P+R*j - k) = XTX(mix);
                    end
                end
            end
            if M > 0
                for j = M-m+1:M
                    %Mixed sums of Laurent base with other polynomial bases.
                    for k = 1:N-n
                        mix = P -R*j + k;
                        XTX(mix) = sum(L_el(1+Ny*(j-1):Ny*j).*P_el(1+Ny*(k-1):Ny*k));
                        XTX(P-j+R*k) = XTX(mix);%Mirror about diagonal
                    end
                    %Mixed sums of Laurent base with other Laurent bases.
                    for k = 1:M-m
                        mix = P-R*j-k;
                        XTX(mix) = sum(L_el(1+Ny*(j-1):Ny*j).*L_el(1+Ny*(k-1):Ny*k));
                        XTX(P-j-R*k) = XTX(mix); %Mirror about diagonal
                    end
                end
            end
        end
    else
        error('Orders of regression are invalid for at least one independent variable.')
    end
    
    % Constraint augmentor
    if Nk > 0 %If there are constraints.
        if i == 1
            XTX(P+sPn+1:P+sPn+Nk) = 1;
            XTX(P+(sPn+1)*R:R:P+(sPn+Nk)*R) = 1;
        end
        if n > 0
            XTX((P+sPn+1:P+sPn+Nk) + R*(N-n+1)) = K(1+Nk*(i-1):Nk*(i));
            XTX((P+R*(sPn+1):R:P+R*(sPn+Nk))+N-n+1) = K(1+Nk*(i-1):Nk*(i));
            for j = N-n+2:N
                col = (sPn+1:sPn+Nk) +P + R*j;
                XTX(col) = XTX(col-R).*K(1+Nk*(i-1):Nk*i);
                XTX((R*(sPn+1):R:R*(sPn+Nk))+ P + j) = XTX(col);
            end
        end
        if m > 0
            XTX((sPn+1:sPn+Nk)+P - R*(M-m+1)) = 1./K(1+Nk*(i-1):Nk*i);
            XTX((R*(sPn+1):R:R*(sPn+Nk))+P - (M-m+1)) = 1./K(1+Nk*(i-1):Nk*i);
            for j = M-m+2:M
                col = (sPn+1:sPn+Nk)+P - R*j;
                XTX(col) = XTX(col+R)./K(1+Nk*(i-1):Nk*i);
                XTX((R*(sPn+1):R:R*(sPn+Nk))+P - j) = XTX(col);
            end
        end
    end
end
coeff = XTX\RHS;
M = 0;
N = 0;
for i = 1:Qx
    p_i = 1 + 2*(i-1);%LI of current variable's polynomial order.
    n = o(p_i);%Current variable's polynomial order.
    m = o(p_i+1);%Current variable's laurent order.
    M = M + m;
    N = N + n;
    OLS.C{i} = [coeff(sLn-M+1:sLn-M+m);coeff(sLn+1+N-n:sLn+N)];
end
data.res = Y - custompolyval(X,OLS.C,o(1,:),coeff(sLn+1));

res2 =Y - (coeff(sLn-3)./(x2.^2)+coeff(sLn-2)./x2 +...
     coeff(sLn-1)./(x1.^2)+coeff(sLn)./x1 +...
     coeff(sLn+1) +...
     coeff(sLn+2)*x1 + coeff(sLn+3)*x1.^2 + coeff(sLn+4)*x1.^3 +...
     coeff(sLn+5)*x3 + coeff(sLn+6)*x3.^2 + coeff(sLn+7)*x3.^3);
 %}
 [data.res,res2]
%% Data Visualization
%{
clc
clear
close all
f = figure('units','normalized');
wfp0 = 0.6;
hfp0 = 0.1;
p0 = uipanel('Parent',f,'Title','Data Visualization',...
    'Position',[1-wfp0,1-hfp0,wfp0,hfp0],...
    'units','normalized');
wfp1 = wfp0;
hfp1 = 1-hfp0;
p1 = uipanel('Parent',f,'Title','Trends',...
    'Position',[1-wfp1,0,wfp1,hfp1],'Visible','on');
p11 = uipanel('Parent',f,'Title','Histograms',...
    'Position',[1-wfp1,0,wfp1,hfp1],'Visible','off');
p12 = uipanel('Parent',f,'Title','Residual Dist.',...
    'Position',[1-wfp1,0,wfp1,hfp1],'Visible','off');
wfp2 = 1-wfp1;
hfp2 = 0.6;
p2 = uipanel('Parent',f,'Title','Descriptive Statistics',...
    'Position',[0,0,wfp2,hfp2],'units','normalized');
p3 = uipanel('Parent',f,'Title','Regression Summary',...
    'Position',[0,hfp2,wfp2,1-hfp2]);
b1 = uicontrol('Parent',p0,'Style','radiobutton',...
    'String','Trends','Units','normalized','Position',[0,0,0.25,1]);
b2 = uicontrol('Parent',p0,'Style','radiobutton',...
    'String','Trends','Units','normalized','Position',[0.25,0,0.25,1]);
b3 = uicontrol('Parent',p0,'Style','radiobutton',...
    'String','Trends','Units','normalized','Position',[0.50,0,0.25,1]);
t1 = uitable('Parent',p2,'Units','normalized','Position',[0,0,1,0.9],...
    'RowName',{'min','Q1','median','Q3','max','std','mean','var.','skew.','kurt.'});
AX = cell(1,Q);
for i = 1:Q
    AX{i} = axes('Parent',f,'Box','on','Boxstyle','full');
end

ceil(sqrt(Q));
switch Q
    case 1 %Univariate regression requires 2D representation.
    case 2 %Bivariate regression requires 3D representation.
    otherwise %Higher order regression requires multiple permutations of 3D visualization.
end
%}