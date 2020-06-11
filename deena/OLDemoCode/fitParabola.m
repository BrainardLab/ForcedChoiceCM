function [error, yPred] = fitParabola(params, xVals, yObs)
% Use with fmincon to find the error associated with a given parabola fit
yPred = params(1) .* (xVals - params(2)).^2 + params(3); 

squaredDiffs = zeros(1,length(yPred)); 
for i = 1:length(yPred)
    squaredDiffs(i) = (yObs(i) - yPred(i))^2;  
end 

error = sqrt(mean(squaredDiffs));
end 