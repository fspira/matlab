function [res,res_tmp] = splineDilateArray(input_img, lookup_image, spline, dilationSize, smoothSize)
    %dilationSize
    %smoothSize
    %nargin
    if nargin < 4
        dilationSize = [2, 5]
    end
    
    if nargin < 5
        smoothSize = 3
    end

    res_tmp = zeros(size(input_img));
    lookup_image_smooth = imfilter(lookup_image, fspecial('gaussian', smoothSize, smoothSize / 4));
    
    for i = 1:length(spline)
       res_tmp(spline(i,2), spline(i,1)) = lookup_image_smooth(spline(i,2), spline(i,1)) ;
    end
    
    res_dil_output = imdilate(res_tmp, strel('disk',dilationSize(1) ,8));
    res_dil_large = imdilate(res_tmp, strel('disk',dilationSize(2),8));
    res = imfilter(res_dil_large, fspecial('gaussian', smoothSize, smoothSize / 4));
    res = res .* (res_dil_output > 0);
end


