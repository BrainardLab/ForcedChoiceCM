halfOnStart = 769 * ones(1, 1024);
halfOnStop = 769 * ones(1, 1024);

ol = OneLight;
%ol.setMirrors(squeeze(halfOnStart)', squeeze(halfOnStop)');
ol.setAll(true);
pause(1);
Snd('Play',sin(0:5000));
ol.setMirrors(squeeze(710 * ones(1, 1024)'),...
    squeeze(710 * ones(1, 1024)'));
pause(1);
Snd('Play',sin(0:5000));
ol.setAll(false); 
