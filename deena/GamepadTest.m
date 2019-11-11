primaryPos = 1;
testPos = 1;
stepModePos = 1; % Start with largest step size
stepModes = [20 5 1];
adjustment_length = 100; 


gamePad = GamePad(); 
fprintf('Starting display loop \n');

while(true)
    key = gamePad.getKeyEvent(); 
    if (~isempty(key))
        switch(key.charCode)
            case 'GP:Y' % Exit program
                break;
            case 'GP:A' % Switch step size mode 
                stepModePos = stepModePos + 1; 
                if stepModePos > length(stepModes)
                    stepModePos = 1; 
                end 
            case 'GP:North' % Scale up test intensity
                testPos = testPos + stepModes(stepModePos);
                if testPos > adjustment_length
                    testPos = adjustment_length;
                end
            case 'GP:South' % Scale down test intensity
                testPos = testPos - stepModes(stepModePos);
                if testPos < 1
                    testPos = 1;
                end
            case 'GP:East' % Move towards p1
                primaryPos = primaryPos + stepModes(stepModePos);
                if primaryPos > adjustment_length
                    primaryPos = adjustment_length;
                end
            case 'GP:West' % Move towards p2
                primaryPos = primaryPos - stepModes(stepModePos);
                if primaryPos < 1
                    primaryPos = 1;
                end
        end
        fprintf('User pressed key. Test intensity = %g, red primary = %g, step size = %g \n',...
            testPos, primaryPos, (stepModes(stepModePos)/100.0));
    end
end 
fprintf('user exited the program'); 