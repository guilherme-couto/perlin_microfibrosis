function [density, total_elements, count_1, count_0] = getPatternDensity(presence)

    % Count how many elements are in the presence, specifying the number of 1 and 0
    total_elements = numel(presence);
    count_1 = sum(presence(:) == 1);
    count_0 = sum(presence(:) == 0);
    density = count_1/total_elements;
    
end