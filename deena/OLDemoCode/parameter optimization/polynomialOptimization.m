function fminconParams = polynomialOptimization(xData, yData, degree, varargin)
% Plots a dataset and finds a polynomial fit using fmincon with standard
% settings. Can enter lower/upper limits and linear equalities/inequalities 
% as key-value pairs. 

% Parse input 
p = inputParser;
p.addParameter('x0',ones(1, degree+1),@(x)(isnumeric(x))); % Initial values 
p.addParameter('A',[],@(x)(isnumeric(x)));    % Linear inequality
p.addParameter('B',[],@(x)(isnumeric(x)));    % Linear inequality
p.addParameter('Aeq',[],@(x)(isnumeric(x)));  % Linear equality
p.addParameter('Beq',[],@(x)(isnumeric(x)));  % Linear equality
p.addParameter('ll',-Inf,@(x)(isnumeric(x))); % Lower limit for all params
p.addParameter('ul',Inf,@(x)(isnumeric(x)));  % Upper limit for all params
p.parse(varargin{:});

% Plot data
figure(); 
hold on; 
plot(xData,yData,'bo '); 

% Lower and upper limits 
lb = p.Results.ll * ones(1,degree+1); 
ub = p.Results.ul * ones(1,degree+1); 

% fmincon settings 
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off',...
    'LargeScale','off','Algorithm','active-set');

% Run optimization 
fminconParams = fmincon(@(x) FitPolynomial(x,xData,yData),...
    p.Results.x0,p.Results.A,p.Results.B,p.Results.Aeq,p.Results.Beq,...
    lb,ub,[],options);

% Generate fit 
yFit = zeros(1,length(yData)); 
for i = 1:(degree+1)
    yFit = yFit + fminconParams(i).*xData.^(degree+1-i); 
end 

% Plot results
plot(xData,yFit,'r-'); 
legend('Data','Fit'); 
theTitle = sprintf('Polynomial Fitting Exercise,degree %g',degree);
title(theTitle);
end 