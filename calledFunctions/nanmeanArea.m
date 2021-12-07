function [out] = nanmeanArea(in,inArea)
% calculates the area weighted mean of pixels without NaN
    inArea = ones(1,size(in,2)) .* inArea;

    idx = find(isnan(in));
    inArea(idx) = NaN;
    
    % total area
    inAreaSum =  sum(inArea, 1, 'omitnan');
    
    % area weighted values
    out = sum(in .* inArea, 1, 'omitnan') ./ inAreaSum;
    
    
end