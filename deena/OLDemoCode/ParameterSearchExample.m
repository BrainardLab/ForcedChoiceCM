% My first attempt to use fmincon to find parameters of a parabola 

%% Generate data: parabola of the form 3(x-1)^2 + 2, with noise added 
xData = -20:1:20; 
yData = 2.5 * ((xData - 3).^2) + 4;
yNoisy = yData + normrnd(0, 4, 1, length(xData)); 

figure; 
hold on; 
plot(xData, yNoisy, 'bo '); 

%% Fit a model of the form y = a(x-b)^2 + c using fmincon
% Initial parameter values and limits 
a = 1; 
b = 1; 
c = 1; 
x0 = [a b c]; 
lb = [-10 -10];
ub = [10 10];

% Options for fmincon 
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off',...
    'LargeScale','off','Algorithm','active-set');
  
% Call optimization function
fminconParams = fmincon(@(x)fitParabola(x,xData,yNoisy),...
    x0,[],[],[],[],lb,ub,[],options);
yFit =  fminconParams(1) * ((xData - fminconParams(2)).^2) + fminconParams(3);
plot(xData, yFit, 'r-'); 
legend('Data', 'Fit'); 
title('Parabola Fitting Exercise'); 