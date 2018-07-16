function [z_scores, bestfocus_z] = bestFocusLevel(imdir, z_stacks, tile, xwell, ywell)
    z_scores = zeros(size(z_stacks));
    wellradius = 110;
    
    files = dir([imdir '*m' sprintf('%03d',tile) '*.tif']);
    
    for i = 1:length(z_stacks)
        z = z_stacks(i);
        try
            IM = double(imread([imdir files(z).name]));
        catch
            continue
        end
        
        [height,width] = size(IM);
        
        % Round up all values <1
        xwell(xwell<1) = 1;
        ywell(ywell<1) = 1;
        xwell(xwell>width) = width;
        ywell(ywell>height) = height;
        
        background = medfilt2(IM, [40 40], 'symmetric');
        foreground = imgaussfilt(IM - background, 2);
        
        % Compute the distance transform.
        BW = zeros(size(IM));
        
        BW(round(ywell), round(xwell)) = 1;
        dist = bwdist(BW);
        
        % Create a mask for all the wells in the image.
        mask = dist < wellradius;
        
        % Mask the foreground image.
        foreground(~mask) = 0;
        
        LOG = fspecial('log', 20, 1);
        ED = imfilter(foreground, LOG, 'symmetric');
        
        thresh = -30;
        crop = imerode(mask, strel('disk', 5));
        BW = crop .* (ED < thresh);
        
%         figure;
%         imshow(ED, [])
%         hold on
%         contour(BW, 'r')
%         title(['stack ' num2str(z)])
        
        whitePixelRatio = sum(BW(:)) / sum(mask(:));
        z_scores(i)= whitePixelRatio;
    end
    
    [~, bestfocus_idx] = max(z_scores);
    bestfocus_z = z_stacks(bestfocus_idx);
end
