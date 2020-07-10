file = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\precomputedStartStops\OLRayleighMatchFinerSpectralSettings_670_560_600_1_0.02_0.07.mat');
primaries = 208:216;
tests = 180:191;
observer = genRayleighObserver();
nCombos = length(primaries)*length(tests);
res = zeros(nCombos,4);
count = 1;
for i = 1:length(primaries)
    for j = 1:length(tests)
        res(count,1) = primaries(i); 
        res(count,2) = tests(j);
        [~,~,isMatch,diffVec] = observerRayleighDecision(observer,...
            file.primarySpdsPredicted(:,primaries(i)),...
            file.testSpdsPredicted(:,tests(j)));
        res(count,3) = isMatch;
        res(count,4) = diffVec; 
    count = count + 1; 
    end 
end 

diffs = res(:,4);
[~,row] = min(diffs); 

%% now, try replicating with single spds - like ideal
res2 = zeros(nCombos,4);
count = 1;
p1Base = makeOLRayleighPrimary(file.p1);
p2Base = makeOLRayleighPrimary(file.p2,'p2',true); 
testBase = makeOLRayleighPrimary(file.test); 
darkSpd = makeOLRayleighPrimary(0);

for i = 1:length(primaries)
    for j = 1:length(tests)
        primaryScale = file.p1Scales(primaries(i));
        testScale = file.testScales(tests(j)); 
        
        testSpd = testBase*testScale*0.07 + darkSpd; 
        primarySpd = p1Base*primaryScale + p2Base*(1-primaryScale)*0.02 + darkSpd;
        
        res2(count,1) = testScale;
        res2(count,2) = primaryScale;
        [~,~,isMatch,diffVec] = observerRayleighDecision(observer,...
            primarySpd,testSpd);
        res2(count,3) = isMatch;
        res2(count,4) = diffVec; 
    count = count + 1; 
    end 
end 

diffs2 = res2(:,4);
[~,row2] = min(diffs2); 