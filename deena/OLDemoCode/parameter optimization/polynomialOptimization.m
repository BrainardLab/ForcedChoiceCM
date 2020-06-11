function fminconParams = polynomialOptimization(xData, yData, degree)
% Plots a dataset and finds a polynomial fit using fmincon

% Plot data
figure(); 
hold on; 
plot(xData, yData, 'bo '); 

% Set all initial parameters to 1 
x0 = ones(1, degree);

% Lower and upper limits 
lb = -10 * x0; 
ub = 10 * x0; 

% fmincon settings 
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off',...
    'LargeScale','off','Algorithm','active-set');

% Run optimization 
[fminconParams, yFit] = fmincon(@(x) FitPolynomial(x,xData,yData),...
    x0,[],[],[],[],lb,ub,[],options);

% Plot results
plot(xData, yFit, 'r-'); 
legend('Data', 'Fit'); 
theTitle = sprintf('Polynomial Fitting Exercise, degree %g', length(fMinConParams)-1);
title(theTitle);
end 