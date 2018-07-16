clear all; close all; clc;

addpath('supp_functions');

tile = 2;
wellnum = 3;
zs = 1:51;

N_stacks = 51;

timepoints = 1:4;
% Loop over timepoints.
for t=timepoints
    % Get the correct image directory name (varies for each timepoint).
    if t == 1
        imdir = 'Y:/Claire/2018-05-17/9717-Image Export-01/';
    else
        imdir = ['Y:/Claire/2018-05-17/9717-' sprintf('%02d', t) ...
            '-Image Export-' sprintf('%02d', t) '/'];
    end
    
    f = dir([imdir '*.tif']);
    filelist={f.name};
    filelist=sort(filelist);
    
    % Loop over all important stacks to find image that will give us the
    % best well centers when passed into coopgerm_wellregister.
    best_focus_value = 0;
    for j = 35:N_stacks
        imfile = char(filelist(tile + 256 * (j-1)));
        % find the tile # and the corresponding z stack for this image
        ind_z = regexp(imfile,'z\d\d');
        is_3digittileid = any(regexp(imfile(ind_z+6),'\d'));
        if is_3digittileid
            z = str2double(imfile(ind_z+(1:2)));
            m = str2double(imfile(ind_z+(4:6)));
        else
            z = str2double(imfile(ind_z+(1:2)));
            m = str2double(imfile(ind_z+(4:5)));
        end
        % Store the best_z_stack for each tile
        im = imread([imdir imfile]);
        im = uint16(im);
        focus_value = focusmeasure(im, 'LAPV');
        % We assume that the computed focus_values is always positive.
        if focus_value >= best_focus_value
            topfocus_im = im;
            topfocus_z = j;
        end
    end
    
    topfocus_z
    
    % Get the well centers.
    [xwell, ywell] = coopgerm_wellregister(topfocus_im);
    
    %figure;
    %imshow(topfocus_im, [])
    %hold on
    %plot(xwell, ywell, 'ro')
    %hold off
    
    xwell = xwell(wellnum);
    ywell = ywell(wellnum);
    
    [z_scores, bestfocus_z] = bestFocusLevel(imdir, zs, tile, xwell, ywell);

    figname = ['m' sprintf('%03d', tile) '_w' sprintf('%02d', wellnum) '_t' sprintf('%02d', t)];
    figure;
    plot(zs, z_scores)
    title(figname)
    savefig(['plots/' figname '.fig'])
end
