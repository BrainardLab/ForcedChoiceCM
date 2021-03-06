% t_AnalyticRayleighMatch
%
% Tutorial/test program to illustrate analytic way of computing Rayleigh
% match and compare with other ways of doing so.
%
% The analytic method does a calculation based on the spectra of the
% red/green primaries and test, and assumes that all mixtures are linear
% combinations of these.
%
% The match can also be computed by exhaustively searching all the
% precomputed mixtures that we use in the OL experimental program. In this
% case, there is quantization error.  The fineness of sampling is
% controlled by parameter adjustmentLength.
%
% When we use the precomputed mixtures, we can generate these either
% nominally, so that they are in fact linear combinations of the nominal
% spectra of the red/green primaries and the test.  Or we can use the
% actual spectra generated by the OL calibration routines, where this
% assumption holds to approximation.  The logical variable NOMINAL controls
% this.
%
% To get all methods to line up, use adjustmentLength of 4001 and set
% NOMINAL to true. This takes a while to run, particularly if you haven't
% already precomputed all the OL spectra. If you shorten the adjustment
% length to 501 you can see small effects of quantization.  If you set
% NOMINAL to true, you can see larger effects of the non-linearty in OL
% spectral generation.
%
% Some care is required to think through and handle dark light, as you need
% to know when it has and hasn't already been added in to the precomputed
% values.  I think it is done right here.
%
% 08/07/2020  dhb  Wrote it.

%% Clear
clear; close all;

%% Define primary and test peak wavelengths
p1Wl = 670;
p2Wl = 560;
testWls = [580 590 600 610 620 630 640 650];
p1ScaleFactor = 1;
p2ScaleFactor = 0.02;
testScaleFactor = 0.5;
adjustmentLength = 501;
FORCE_COMPUTE_OLLIGHT = false;
NOMINAL = false;

%% Loop over test peak wavelengths
for ww = 1:length(testWls)
    % Get test wavelength
    testWl = testWls(ww);
    
    % Get OL spectra, precmputing if needed or forced
    startsStopsFilename = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
        adjustmentLength, p1Wl, p2Wl, testWl, p1ScaleFactor, p2ScaleFactor, testScaleFactor);
    if (FORCE_COMPUTE_OLLIGHT || ~exist(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'precomputedStartStops', startsStopsFilename),'file'))
        fprintf('Precomputing starts/stops\n');
        OLRayleighMatchLightSettings(p1Wl, p2Wl, testWl, ...
            'p1ScaleFactor',p1ScaleFactor,'p2ScaleFactor',p2ScaleFactor,'testScaleFactor',testScaleFactor, ...
            'adjustmentLength',adjustmentLength);
    end
    
    % Load the various OL spectra and related info.
    fprintf('Loading precomputed starts/stops\n');
    olData = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'precomputedStartStops', startsStopsFilename));
    S = WlsToS(olData.wls);
    
    % Some plots of spectra
    %
    %   olData.primarySpdsPredicted - Predicted mixtures of primaries available in the experiment,includes dark light.
    %   olData.primary1IncrSpd      - Nominal primary 1 at full power, dark light subtracted
    %   olData.primary2IncrSpd      - Nominal primary 2 at full power, dark light subtracted
    %   olData.primary1IncrSpdPredicted - Predicted primary 1 at full power, dark light subtracted
    %   olData.primary2IncrSpdPredicted - Predicted  primary 2 at full power, dark light subtracted
    if (~exist('primaryFig','var'))
        primaryFig = figure; clf; hold on
    else
        figure(primaryFig); clf; hold on
    end
    for ii = 1:size(olData.primarySpdsPredicted,2)
        plot(olData.wls,olData.primarySpdsPredicted(:,ii)-olData.darkSpd,'c','LineWidth',0.5);
    end
    plot(olData.wls,p1ScaleFactor*olData.primary1IncrSpdPredicted,'r','LineWidth',3);
    plot(olData.wls,p2ScaleFactor*olData.primary2IncrSpdPredicted,'g','LineWidth',3);
    plot(olData.wls,p1ScaleFactor*olData.primary1IncrSpd,'k:','LineWidth',2);
    plot(olData.wls,p2ScaleFactor*olData.primary2IncrSpd,'k:','LineWidth',2);

    %    olData.testSpdsPredicted  - Predicted test lights available in the experiment, includes dark light.
    %    olData.testIncrSpd        - Nominal test light at full power, dark light subtracted
    %    olData.testIncrSpd        - Predicted test light at full power, dark light subtracted
    %    olData.testIncrSpdPredicted - Predicted test light at full power, dark light subtracted
    if (~exist('testFig','var'))
        testFig = figure; clf; hold on
    else
        figure(testFig); clf; hold on
    end
    for ii = 1:size(olData.testSpdsPredicted,2)
        plot(olData.wls,olData.testSpdsPredicted(:,ii)-olData.darkSpd,'c','LineWidth',0.5);
    end
    plot(olData.wls,testScaleFactor*olData.testIncrSpdPredicted,'y','LineWidth',3);
    plot(olData.wls,testScaleFactor*olData.testIncrSpd,'k:','LineWidth',3);
    
    % Set up whether we use nominal or predicted in everything below.
    % Using the NOMINAL values leads to performance consistent with theory,
    % as long as the value of adjustmentLength is large.
    if (NOMINAL)
        testSpds = olData.testSpdsNominal;
        primarySpds = olData.primarySpdsNominal;
        testIncrSpd = olData.testIncrSpd;
        primary1IncrSpd = olData.primary1IncrSpd;
        primary2IncrSpd = olData.primary2IncrSpd;
    else
        testSpds = olData.testSpdsPredicted;
        primarySpds = olData.primarySpdsPredicted;
        testIncrSpd = olData.testIncrSpdPredicted;
        primary1IncrSpd = olData.primary1IncrSpdPredicted;
        primary2IncrSpd = olData.primary2IncrSpdPredicted;
    end
    
    % Get an observer
    observer = genRayleighObserver('S',S);
    observerParamsVec = ObserverParamsToVec('basic',observer);
    opponentParamsVec = [observer.colorDiffParams.lumWeight, observer.colorDiffParams.rgWeight, ...
        observer.colorDiffParams.byWeight, observer.colorDiffParams.noiseSd];
    T_cones = observer.T_cones;
    
    % Compute analytic predicted Rayleigh match. Use predicted incremtental
    % stimuli so we can compare with stimuli available in experiment below.
    % Also pass zero dark light to keep computedPredictedRayleigh() from
    % doing anything fancy and confusing us.  We handle dark light
    % addition and subtraction explicitly at this level.
    [testAdjustSpd,primaryMixtureSpd,testIntensity(ww),lambda(ww)] = ...
        computePredictedRayleighMatch(p1ScaleFactor*primary1IncrSpd, ...
        p2ScaleFactor*primary2IncrSpd, ...
        testScaleFactor*testIncrSpd, ...
        observerParamsVec(1:8),opponentParamsVec, ...
        'fieldSize',observer.coneParams.fieldSizeDegrees,'age',observer.coneParams.ageYears,'S',S, ...
        'darkSpd',zeros(size(olData.darkSpd)));
    testMatchSpd = testAdjustSpd + olData.darkSpd;
    primaryMatchSpd = primaryMixtureSpd + olData.darkSpd;
    
    % Check that LM cone excitations match to high precision
    testMatchLM(:,ww) = T_cones(1:2,:)*testMatchSpd;
    primaryMatchLM(:,ww) = T_cones(1:2,:)*primaryMatchSpd;
    if (max(abs( (testMatchLM(:,ww)-primaryMatchLM(:,ww))./testMatchLM(:,ww) )) > 1e-6)
        error('Predicted Rayleigh match does not equate LM cone exciations');
    end
    
    % Add Rayleigh match spectra to plots
    figure(testFig);
    plot(olData.wls,testAdjustSpd,'k','LineWidth',3);
    figure(primaryFig);
    plot(olData.wls,primaryMixtureSpd,'k','LineWidth',3);
    
    % Find closest match to predicted spectra in availble tests and primary
    % mixtures. Just do exhaustive search to find closest in each table.
    % This is just as a check, not all that useful in the end.  More useful
    % and more diagnostic is the exhaustive search below.
    minDiff = Inf;
    index = NaN;
    for tt = 1:size(testSpds,2)
        diff = sum((testMatchSpd-testSpds(:,tt)).^2);
        if (diff < minDiff)
            minDiff = diff;
            index = tt;
        end
    end
    testIndex(ww) = index;
    searchTestIntensity(ww) = testIndex(ww) / size(testSpds,2);
    testClosestSpd = testSpds(:,testIndex(ww));
    
    minDiff = Inf;
    index = NaN;
    for pp = 1:size(primarySpds,2)
        diff = sum((primaryMatchSpd-primarySpds(:,pp)).^2);
        if (diff < minDiff)
            minDiff = diff;
            index = pp;
        end
    end
    primaryIndex(ww) = index;
    searchLambda(ww) = primaryIndex(ww) / size(primarySpds,2);
    primaryClosestSpd = primarySpds(:,primaryIndex(ww));
    
    %% Check that LM cone excitations match well enough for closest approximation
    %
    % Whether they do or not depends on both NOMINAL and the number of
    % adjustmentLevels, so the actual error check is commented out.
    testClosestLM(:,ww) = T_cones(1:2,:)*testClosestSpd;
    primaryClosestLM(:,ww) = T_cones(1:2,:)*primaryClosestSpd;
    if (max(abs( (testClosestLM(:,ww)-primaryClosestLM(:,ww))./testClosestLM(:,ww) )) > 1e-2)
        %error('Predicted closest Rayleigh match does not equate LM cone exciations to 1%%');
    end
    
    %% Add closest approximation to plots
    figure(testFig);
    plot(olData.wls,testClosestSpd-olData.darkSpd,'r:','LineWidth',2);
    figure(primaryFig);
    plot(olData.wls,primaryClosestSpd-olData.darkSpd,'r:','LineWidth',2);
    
    % Do exhaustive search to find available test/primary pair that
    % minimizes LM difference.  This might usefully be made a funciton to
    % go with the OL simulation suite, as it is the method that always
    % tells us the ideal Rayleigh match when the spectra are quantized and
    % linearity does not hold perfectly.  Error check is commented out
    % because it only works for some choices of NOMINAL and
    % adjustmentLevels.
    minDiff = Inf;
    ttIndex(ww) = NaN;
    ppIndex(ww) = NaN;
    LMSTestTemp = T_cones*testSpds;
    LMSPrimaryTemp = T_cones*primarySpds;
    for tt = 1:size(testSpds,2)
        opponentDiffTemp = LMSToOpponentContrast(observer.colorDiffParams,LMSTestTemp(:,tt),LMSPrimaryTemp);
        diffTemp = vecnorm(opponentDiffTemp(1:2,:));
        [diff,pp] = min(diffTemp);
        if (diff < minDiff)
            minDiff = diff;
            ttIndex(ww) = tt;
            ppIndex(ww) = pp;
            ttLMTest(:,ww) = LMSTestTemp(1:2,tt);
            ppLMPrimary(:,ww) = LMSPrimaryTemp(1:2,pp);
            ttSpdTest(:,ww) = testSpds(:,tt);
            ppSpdPrimary(:,ww) = primarySpds(:,pp);
        end
    end
    ttTestIntensity(ww) = ttIndex(ww) / size(testSpds,2);
    ppLambda(ww) = ppIndex(ww) / size(primarySpds,2);
    if (abs(ttIndex(ww)-testIndex(ww)) > 1)
        %error('Don''t recover same test index within rouding in two ways');
    end
    if (abs(ppIndex(ww)-primaryIndex(ww)) > 1)
        %error('Don''t recover same primary index witin rounding in two ways');
    end   
    
    % Check that analytic indices are within rounding.  Another one where
    % the error call is commented out.
    analyticTestIndex(ww) = round(size(testSpds,2)*testIntensity(ww));
    analyticPrimaryIndex(ww) = round(size(primarySpds,2)*lambda(ww));
    if (abs(analyticTestIndex(ww)-ttIndex(ww)) > 1)
        %error('Analytic test index not within rouding of exhaustive search');
    end
    if (abs(analyticPrimaryIndex(ww)-ppIndex(ww)) > 1)
        %error('Analytic primary index not within rouding of exhaustive search');
    end
    
    % Check whether primaries are linear combination of primary1 and 2 and
    % whether test is a scalar multiply of the test.  Close but not
    % exactly. It turns out that the small differences are quite salient
    % when we plot their effect in the generalized Pitt diagram.
    ttSpdIncrTestFit = testIncrSpd*(testIncrSpd\(ttSpdTest(:,ww)-olData.darkSpd));
    ppSpdIncrPrimaryFit = [primary1IncrSpd primary2IncrSpd]*([primary1IncrSpd primary2IncrSpd]\(ppSpdPrimary(:,ww)-olData.darkSpd));
    figure; clf;
    subplot(2,2,1); hold on
    plot(olData.wls,ttSpdTest(:,ww)-olData.darkSpd,'k','LineWidth',5);
    plot(olData.wls,ttSpdIncrTestFit,'r','LineWidth',3);
    subplot(2,2,2); hold on
    plot(olData.wls,ppSpdPrimary(:,ww)-olData.darkSpd,'k','LineWidth',5);
    plot(olData.wls,ppSpdIncrPrimaryFit,'g','LineWidth',3);
    subplot(2,2,3); hold on
    plot(ttSpdTest(:,ww)-olData.darkSpd,ttSpdIncrTestFit,'r+','MarkerSize',8);
    plot([0 max([ttSpdTest(:,ww)-olData.darkSpd ; ttSpdIncrTestFit])],[0 max([ttSpdTest(:,ww)-olData.darkSpd ; ttSpdIncrTestFit])],'k','LineWidth',0.5);
    xlim([0 max([ttSpdTest(:,ww)-olData.darkSpd ; ttSpdIncrTestFit])]);
    ylim([0 max([ttSpdTest(:,ww)-olData.darkSpd ; ttSpdIncrTestFit])]);
    subplot(2,2,4); hold on
    plot(ppSpdPrimary(:,ww)-olData.darkSpd,ppSpdIncrPrimaryFit,'g+','MarkerSize',8);
    plot([0 max([ppSpdPrimary(:,ww)-olData.darkSpd ; ppSpdIncrPrimaryFit])],[0 max([ppSpdPrimary(:,ww)-olData.darkSpd ; ppSpdIncrPrimaryFit])],'k','LineWidth',0.5);
    xlim([0 max([ppSpdPrimary(:,ww)-olData.darkSpd ; ppSpdIncrPrimaryFit])]);
    ylim([0 max([ppSpdPrimary(:,ww)-olData.darkSpd ; ppSpdIncrPrimaryFit])]);
   
end

%% Plot generalized Pitt diagram
pittDiagram = figure; clf; hold on;
plot(lambda,testIntensity,'ro','MarkerFaceColor','r','MarkerSize',20);
plot(searchLambda,searchTestIntensity,'ko','MarkerFaceColor','k','MarkerSize',18);
plot(ppLambda,ttTestIntensity,'ko','MarkerFaceColor','g','MarkerSize',16);