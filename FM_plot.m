wellfile = 'Y:\Hl_wells\2018-03-07_wells\m001-well08.tif';

    roi = false(221);
    roi(111,111) = 1;
    roi = imdilate(roi,strel('disk',37,0));
    FM = zeros(6,1);
    figure;
for i_z = 1:6
    im = imread(wellfile,i_z);
    if i_z==5
        imshow(imfuse(im,edge(im)));
    end
    FM(i_z) = focusmeasure(im,roi);
end

plot(FM, 'r');
set(gca, 'XTick', 1:6)
xlabel('#pic');
ylabel('variance');
title('Focus Measure');

function FM = focusmeasure(Image, ROI)
    %This function measures the relative degree of focus of
    %an image. It may be invoked as:
    %
    %   FM = focusmeasure(IMAGE, ROI)
    %
    %Where
    %   IMAGE,  is a grayscale image and FM is the computed
    %           focus value.
    %   ROI,    Image ROI as a rectangle [xo yo width heigth].
    %           if an empty argument is passed, the whole
    %           image is processed.

    if nargin>2 && ~isempty(ROI)
        Image = imcrop(Image, ROI);
    end
    LAP = fspecial('laplacian');
    ILAP = imfilter(Image, LAP, 'replicate', 'conv');
    FM = std2(ILAP)^2;
end

    