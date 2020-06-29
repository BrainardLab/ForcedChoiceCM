function error = FitPolynomial(params, xVals,yObs)
% Use with fmincon to find the error associated with a given polynomial fit

% Find the predicted y values
yPred = zeros(1, length(yObs)); 
degree = length(params); 
for i = 1:degree
    yPred = yPred + params(i).*xVals.^(degree-i); 
end

% Calculate the root mean square error
squaredDiffs = zeros(1,length(yPred));
for j = 1:length(yPred)
    squaredDiffs(j) = (yObs(j) - yPred(j))^2;  
end 
error = sqrt(mean(squaredDiffs));
end 