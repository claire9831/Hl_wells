%% Focus measure using the variance of laplacian method
function FM = focusmeasure(Image, ROI)
    %This function measures the relative degree of focus of
    %an image.
    if nargin>2 && ~isempty(ROI)
        Image = imcrop(Image, ROI);
    end
    LAP = fspecial('laplacian');
    ILAP = imfilter(Image, LAP, 'replicate', 'conv');
    FM = std2(ILAP)^2;
end
