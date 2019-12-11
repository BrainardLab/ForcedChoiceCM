halfOnStart = zeros(1, 1024);
halfOnStop = zeros(1, 1024);
for i = 1:1024
    halfOnStart(i) = 384;
    halfOnStop(i) = 384;
end
ol = OneLight;
ol.setAll(false);
pause(1);
ol.setMirrors(squeeze(halfOnStart)', squeeze(halfOnStop)');
pause(1);
ol.setAll(true);