%% Find the best match among a set of spectra
% Use OL spectra
close all;
file = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\precomputedStartStops\OLRayleighMatch1601SpectralSettings_670_560_600_1_0.02_0.5.mat');

observer = genRayleighObserver('coneVec',[0 0 0 0 0 0 0 0 0]);

[~,nSpds] = size(file.primarySpdsNominal);
inds = 1:nSpds;
nCombos = nSpds^2;
res = zeros(nCombos,4);
count = 1;

for i = 1:length(inds)      % Loop through primary mixtures
    for j = 1:length(inds)  % Loop through test lights
        res(count,1) = inds(i);
        res(count,2) = inds(j);
        primarySpd = file.primarySpdsNominal(:,inds(i));
        testSpd = file.testSpdsNominal(:,inds(j));
        [~,~,isMatch,diffVec] = observerRayleighDecision(observer,...
            primarySpd,testSpd,'thresholdScale',0.5);
        res(count,3) = isMatch;
        res(count,4) = diffVec;
        count = count + 1;
    end
end

% Find best match
diffs = res(:,4);
[~,row] = min(diffs);

% Plot all matches
matches = res((res(:,3)==1),:);
[nMatches,~] = size(matches);
figure(1);
title('Matches Below Threshold');
xlim([0 1]);
ylim([0 1]);
hold on;
for i = 1:nMatches
    plot(file.p1Scales(matches(i,1)),file.testScales(matches(i,2)),'gs ',...
        'MarkerFaceColor','g');
    %         if (mod(i,5)==0)
    %             primaryLMS = observer.T_cones*file.primarySpdsNominal(:,matches(i,1));
    %             testLMS = observer.T_cones*file.testSpdsNominal(:,matches(i,2));
    %             OLPlotConeEffects(primaryLMS,testLMS,'',1);
    %         end
end
[idealTest,idealPrimary,idealTI,idealLambda] = ...
    computePredictedRayleighMatch(670,560,600,zeros(1,9),'p1Scale',1,...
    'p2Scale',0.02,'testScale',0.5);
plot(idealLambda,idealTI,'b* ');

% Plot cone responses of nominal match and of best predicted match
[~,~,isIdealMatch,idealDiffVec] = observerRayleighDecision(observer,...
    file.primarySpdsNominal(:,res(row,1)),file.testSpdsNominal(:,res(row,2)));
primaryLMS = observer.T_cones*file.primarySpdsNominal(:,res(row,1));
testLMS = observer.T_cones*file.primarySpdsNominal(:,res(row,1));
OLPlotConeEffects(primaryLMS,testLMS,'',1);
title('Best Opponent Contrast Match');

[~,~,isIdealMatch2,idealDiffVec2] = observerRayleighDecision(observer,...
    idealPrimary,idealTest);
primaryLMS2 = observer.T_cones*idealPrimary;
testLMS2 = observer.T_cones*idealTest;
OLPlotConeEffects(primaryLMS2,testLMS2,'',1);
title('Computed Cone Response Match');

[~,pIdealIndex] = min(abs(file.p1Scales-idealLambda));
[~,tIdealIndex] = min(abs(file.testScales-idealTI));
[~,~,isIdealMatch3,idealDiffVec3] = observerRayleighDecision(observer,...
    file.primarySpdsNominal(:,pIdealIndex),...
    file.testSpdsNominal(:,tIdealIndex),'thresholdScale',0.5);
primaryLMS3 = observer.T_cones*idealPrimary;
testLMS3 = observer.T_cones*idealTest;
OLPlotConeEffects(primaryLMS3,testLMS3,'',1);
title('Closest to Ideal Match Cone Responses');
%% now, try replicating with single spds - like ideal
res2 = zeros(nCombos,4);
count = 1;
p1Base = makeOLRayleighPrimary(file.p1);
p2Base = makeOLRayleighPrimary(file.p2,'p2',true);
testBase = makeOLRayleighPrimary(file.test);
darkSpd = makeOLRayleighPrimary(0);

for i = 1:length(primaryInds)
    for j = 1:length(testInds)
        primaryScale = file.p1Scales(primaryInds(i));
        testScale = file.testScales(testInds(j));
        
        testSpd = testBase*testScale*0.5 + darkSpd;
        primarySpd = p1Base*primaryScale + p2Base*(1-primaryScale)*0.02 + darkSpd;
        
        res2(count,1) = testScale;
        res2(count,2) = primaryScale;
        [~,~,isMatch,diffVec] = observerRayleighDecision(observer,...
            primarySpd,testSpd,'thresholdScale',0.5);
        res2(count,3) = isMatch;
        res2(count,4) = diffVec;
        count = count + 1;
    end
end

diffs2 = res2(:,4);
[~,row2] = min(diffs2);

%% Spd scaling tricks
file = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\precomputedStartStops\OLRayleighMatchFinerSpectralSettings_670_560_600_1_0.02_0.5.mat');
% Base spds 
p1Base = makeOLRayleighPrimary(file.p1);
p2Base = makeOLRayleighPrimary(file.p2,'p2',true);
testBase = makeOLRayleighPrimary(file.test);
darkSpd = makeOLRayleighPrimary(0);

% Convert to predicted at the end 
primaryScale = 0.75;
testScale = 0.25;
testSpdNominal = testBase*testScale + darkSpd;
primarySpdNominal = p1Base*primaryScale + p2Base*(1-primaryScale)...
    + darkSpd;
[~,~,primarySpdScaleEnd] = OLSpdToSettings(file.cal,primarySpdNominal,...
    'lambda',file.lambda);
[~,~,testSpdScaleEnd] = OLSpdToSettings(file.cal,testSpdNominal,...
    'lambda',file.lambda);

% Convert to predicted at the beginning 
[~,~,p1BaseScale] = OLSpdToSettings(file.cal,p1Base,'lambda',file.lambda);
[~,~,p2BaseScale] = OLSpdToSettings(file.cal,p2Base,'lambda',file.lambda);
[~,~,testBaseScale] = OLSpdToSettings(file.cal,testBase,'lambda',file.lambda);

primarySpdScaleEarly = p1BaseScale*primaryScale + p2BaseScale...
    *(1-primaryScale) + darkSpd;
testSpdScaleEarly = testBaseScale*testScale + darkSpd;

% Compare the two methods 
OLPlotSpdCheck(SToWls([380 2 201]),[testSpdScaleEnd testSpdScaleEarly]);
title('Test spds');
legend('Scale End','Scale Early'); 

OLPlotSpdCheck(SToWls([380 2 201]),[primarySpdScaleEnd primarySpdScaleEarly]);
title('Primary spds');
legend('Scale End','Scale Early');