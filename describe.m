%% Descriptive Statistics
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE), College of Arts and Sciences (COAS)
%  For use in MA 412, MA 413, AE 435, AE 440 and any other course that
%  would benefit from a data-fitting tool.
%% Description
% Computes the summary statistics of data sets. Three types of statistics
% are considered: Those that measure spread (variance, standard deviation,
% range, those that measure location (mean and median), and those that
% measure shape (skewness and kurtosis).
%% Formulae
%  Arithmetic mean
%%
% $\bar{X} = \frac{1}{N}\left(\sum_{i=1}^{N}{x_{i}}\right)$
%%
%  Variance
%% 
% $Var(x) = \frac{1}{N}\left(\sum_{i=1}^{N}(x_{i} - \bar{X})^{2}\right)$
%% Required Plugins
% * none
%% Changelog
%  v1.0,(08/07/2022): Initial Release. Reinvented the wheel once again.
%% Syntax
% * INPUT(*X*): A numeric array of arbitrary size that conveys the
% observations. The columns are taken to be different variables whereas the
% rows are taken to be the observations.
% * OUTPUT(*summary*): A structure with fields for the  descriptive
% statistics.
%% Function definition
function summary = describe(X)
[rx,cx] = size(X);
summary = struct(...
    'N',ones(1,cx)*rc,... %Observations
    'min',zeros(1,cx),...%Minimum values (zeroth quartile).
    'Q1',zeros(1,cx),... %First quartile.
    'med',zeros(1,cx),...%Median (second quartile).
    'Q3',zeros(1,cx),... %Third quartile.
    'max',zeros(1,cx),...%Maximum values (fourth quartile).
    'ran',zeros(1,cx),...%Range.
    'avg',zeros(1,cx),...%Arithmetic means.
    'std',zeros(1,cx),...%Standard deviation.
    'var',zeros(1,cx),...%Variance.
    'skw',zeros(1,cx),...%Skewness.
    'kur',zeros(1,cx));  %Kurtosis.
for i = 1:cx
    X(:,i) = sort(X(:,i)); %Replace later with a custom script.
    summary.min(i) = X(1,i);%0th-quartile.
    summary.max(i) = X(rx,i);%100th-quartile.
    summary.ran(i) = X(rx,i) - X(1,i);%Range.
    summary.avg(i) = sum(X(:,i))/rx;%Arithmetic mean.
    summary.var(i) = sum((X(:,i)-summary.avg(i)).^2)/rx;%Variance.
    summary.std(i) = sqrt(summary.var(i));%Standard deviation.
    summary.skw(i) = sum((X(:,i)-summary.avg(i)).^3)/(rx*summary.var(i)^(3/2)); %Skewness.
    summary.kur(i) = sum((X(:,i)-summary.avg(i)).^4)/(rx*summary.var(i)^2);%Kurtosis.
    
    %50th-percentile (median).
    if mod(rx,2) == 0 %Even number of observations.
        half = rx/2;
        summary.med(i) = (X(half,i) + X(half + 1,i))/2;
    else %Odd number of observations.
        half = rx+1/2;
        summary.med(i) = X(half,i);
    end
    
    %25th and 75th percentiles (quartiles)
    if mod(half,2) == 0
        summary.Q1(i) = X(half/2,i);
        summary.Q3(i) = X(half + half/2,i);
    else
        summary.Q1(i) = X((half+1)/2,i);
        summary.Q3(i) = X((rx-half + 1)/2,i);
    end
        
end
end