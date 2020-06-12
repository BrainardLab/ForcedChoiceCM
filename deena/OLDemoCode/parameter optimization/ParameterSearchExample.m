% Examples of using the polynomialOptimization function to fit a dataset 
% with polynomials of various degrees 

xData = -20:1:20; % General x data to use throughout the tutorial 

%% Case 1 - Parabola 
% (4/5)x^2 - (2/3)x + 2
% No-noise case 
y1 = (4/5).*xData.^2 - (2/3).*xData + 2; 
parabolaParams = polynomialOptimization(xData,y1,2);
errParabola = FitPolynomial(parabolaParams,xData,y1); 

% With noise added - error is much larger
y1Noisy = y1 + normrnd(0,1,1,length(y1));
parabolaParamsNoisy = polynomialOptimization(xData,y1Noisy,2); 
errParabolaNoisy = FitPolynomial(parabolaParamsNoisy,xData,y1Noisy);

%% Case 2 - Quartic Function 
% 2.5x^4 + 3.2x^3 - 0.2x^2 + 10x - 7.1
% Normal case 
y2 = 2.5*xData.^4 + 3.2*xData.^3 - 0.2*xData.^2 + 10*xData - 7.1; 
quarticParams = polynomialOptimization(xData,y2,4);
errQuartic = FitPolynomial(quarticParams,xData,y2);

% Misfit to wrong degree - too high 
% extra param is effectively set to 0, but error is still higher 
quarticParams2 = polynomialOptimization(xData,y2,5);
errQuarticHigh = FitPolynomial(quarticParams2,xData,y2);

% Misfit to wrong degree - too low 
% Fit is really bad when params are constrained -10 to 10, but better when
% they are unconstrained (still much worse than the other two) 
quarticParams3 = polynomialOptimization(xData,y2,3);
errQuarticLow = FitPolynomial(quarticParams3,xData,y2);

%% Case 3 - Several fits for one dataset 
% Dataset with small x^5 term - can set with 5 or 6 params 
% Fit looks good for degree 4, but the error is much higher
y3 = 0.0025*xData.^5 + 2*xData.^4 - 0.2*xData.^3 - xData.^2 +...
    (2/3)*xData - 0.14; 

degree4Params = polynomialOptimization(xData,y3,4);
errDegree4 = FitPolynomial(degree4Params,xData,y3);

degree5Params = polynomialOptimization(xData,y3,5);
errDegree5 = FitPolynomial(degree5Params,xData,y3);

%% Case 4 - Upper/lower bounds and their violation 
% (4/5)x^2 - 15x + 2
y4 = (4/5)*xData.^2 - 15*xData + 2; 

% Plot with lower limits set to -10 (not low enough)
% Fit is a lot worse
y4Params1 = polynomialOptimization(xData,y4,2,'ll',-10);
errDegree41 = FitPolynomial(y4Params1,xData,y4);

% Plot with lower limits set to -20 (low enough) 
y4Params2 = polynomialOptimization(xData,y4,2,'ll',-20);
errDegree42 = FitPolynomial(y4Params2,xData,y4);

% Plot without constraints 
% Finds very similar params as above, but error is slightly higher (why?)
y4Params3 = polynomialOptimization(xData,y4,2,'ll',-Inf);
errDegree43 = FitPolynomial(y4Params3,xData,y4);

%% Cases 5-6 - Mess up intial values 
% In these cases, the program found the true parameters. However, need to
% be careful for cases where there are local minima
% Try being several orders of magnitude too high 
y5 =  0.8*xData.^2 - 0.5*xData + 0.34; 
y5ParamsHigh = polynomialOptimization(xData,y5,2,'x0',[1000 1000 1000]);
errHigh = FitPolynomial(y5ParamsHigh,xData,y5); 

% Try being several orders of magnitude too low 
y6 =  800*xData.^2 - 0.5*xData + 1000;
y6ParamsLow = polynomialOptimization(xData,y6,2);
errLow = FitPolynomial(y6ParamsLow,xData,y6); 
%% Cases 7-8 - Linear equality/inequality constraints
y7 = 0.5*xData.^3 - 5*xData.^2 + 2*xData.^2 + 0.7;  

% Case 1: all params must sum to 1 
% (this did not find a good fit for this dataset)
y7ParamsEq = polynomialOptimization(xData,y7,3,'Aeq',[1 1 1 1],'Beq',1);
errEq = FitPolynomial(y7ParamsEq,xData,y7);

% Case 2: all params must sum to <= 2 
y7ParamsIneq = polynomialOptimization(xData,y7,3,'A',[1 1 1 1],'B',2);
errIneq = FitPolynomial(y7ParamsIneq,xData,y7);