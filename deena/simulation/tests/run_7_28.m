sampleRayleighMatchSeries('varyNMatchesNoisy3',20,670,560,570:5:640,...
    [0 0 1 1 0 1 1 0],'nMatches',[1 2 3 4],'baseConeParams',...
    [0 0 0 0 0 0 0 0 0.06]);
sampleRayleighMatchSeries('varyIncNMatchesNoisy',20,670,560,570:5:640,...
    [0 0 1 1 0 1 1 0],'testWlIncr',[40 20 10 5 2],...
    'testingParamToVary2','nMatches','testingValsToVary2',[1 2 3],...
    'baseConeParams',[0 0 0 0 0 0 0 0 0.02]);
