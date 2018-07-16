% Created bu Rui Yan, Jan 2018

function Hlwell_tilefocusNew(exper_date, N_tiles, N_stacks, mlist)
% This code finds the best focused images for each well per tile.
% Move to the directory of experimental images

if nargin == 3
    mlist = 1:N_tiles;
end

timepoint = 1;
imdir = [exper_date '/9717-Image Export-0' num2str(timepoint) '/'];

% Get a list of image files in imdirs
f = dir([imdir '*.tif']);
filelist={f.name};
filelist=sort(filelist);


%% Find the best-focused images (best z_stack) for the top of each tile
topfocus_index = ones(N_tiles,1);
topfocus_value = zeros(N_tiles,1);
topfocus_z = ones(N_tiles,1);

% Only loop through the images with z_stack from 35 to 51
for i = mlist
    for j = 35:N_stacks
        imfile = char(filelist(((i-1)*N_stacks)+j));
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
        if focus_value >= best_focus_value(m)
            topfocus_index(m) = (i-1)*N_stacks+j;
            topfocus_value(m) = focus_value;
            topfocus_z(m) = z;
        end
    end
end


%% Get the wells image (contains spores (top)) on each tile
% Also detect the existence of spores on each well image (determine growth)
exptdir = [exper_date, '_wells_t', int2str(timepoint), '/']; 

if 7 ~= exist(exptdir,'dir')
    mkdir(exptdir);
end


% Crop the well image of the top of each well that contains spores
wellradius = 110;
for m = mlist % or list of tiles that we wish to analyze
    topfocus_imfile = char(filelist(topfocus_index(m)));
    topfocus_im = imread([imdir topfocus_imfile]);
    [height, width] = size(topfocus_im);
    [xwell, ywell] = (coopgerm_wellregister(topfocus_im));
    
    bestfocus_z = bestFocusLevel(imdir, timepoint, 15:35, m, xwell, ywell);
    numOfWells = length(xwell);
    scores = cell(numOfWells, 3);
    
    % Read in the tile image.
    imfile = [imdir '9717-Image Export-0' num2str(timepoint) '_z' ...
            num2str(bestfocus_z,'%02i') 'm' num2str(m,'%03i') '_ORG.tif'];
    im = imread(imfile);
    
    
    for i_well = 1:numOfWells
        % Crop each well.
        xrange = round(xwell(i_well))+(-wellradius:wellradius);
        xrange(xrange<1 | xrange>width) = [];
        yrange = round(ywell(i_well))+(-wellradius:wellradius);
        yrange(yrange<1 | yrange>height) = [];
        im_well = im(yrange,xrange);
        
        
        % Save the well image.
        well_name = ['well_' sprintf('%02d', timepoint) ...
                '_z' sprintf('%02d', i) 'm' sprintf('%03d',tile) ...
                'wi' sprintf('%02d', i_well)];
        im_wellfilename = [exptdir well_name '.tif'];
        imwrite(im_well, im_wellfilename);
    end
end

end
