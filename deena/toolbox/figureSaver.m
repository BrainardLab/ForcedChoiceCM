function figureSaver
h =  findobj('type','figure');
n = length(h);
for i = 1:n
    figure(i);
    name = 'c3_od_split.tiff'; 
    fName = fullfile('C:', 'Users', 'deena', 'Documents', 'brainard lab code',...
        'ForcedChoiceCM', 'deena', 'Pitt diagram tutorials', 'pittPlots',...
        name); 
    saveas(gcf, fName); 
end 
end