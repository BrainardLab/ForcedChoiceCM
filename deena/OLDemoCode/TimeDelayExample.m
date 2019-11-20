%% Illustrate how to time things

%% Get current time
delaySeconds = 1;
secondsPerLoop = 3;

%% Little loop that waits
fprintf('Waiting %d seconds at a time in a loop.\n',delaySecondsmg);
nowTime = mglGetSecs;
for ii = 1:secondsPerLoop
    while (mglGetSecs < nowTime + delaySeconds)
    end
    fprintf('\tJust waited %d seconds\n',ii*delaySeconds);
    nowTime = nowTime + delaySeconds;
end
fprintf('Done with first method.\n');

fprintf('Doing it again.\n');
for ii = 1:secondsPerLoop
    mglWaitSecs(delaySeconds);
    fprintf('\tJust waited %d seconds\n',ii*delaySeconds);
end
fprintf('Done with second method.\n');
