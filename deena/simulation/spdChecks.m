calcSpdFile = load('dat1.mat');
[spdLength,nMatches] = size(calcSpdFile.p);
wls = SToWls([380 2 201]); 
for i = 1:nMatches
    file = sprintf('idealTesting_%g.mat',i);
    fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'idealTesting',file);
    matchData = load(fName); 
    
    OLPlotSpdCheck(wls,[matchData.idealPrimarySpd,calcSpdFile.p(:,i)]);
    title1 = sprintf('Primaries,test = %g',matchData.test);
    title(title1);
    legend('Ideal','Simulated');
    
    OLPlotSpdCheck(wls,[matchData.idealTestSpd,calcSpdFile.t(:,i)]);
    title2 = sprintf('Test,test = %g',matchData.test);
    title(title2);
    legend('Ideal','Simulated');
end 
