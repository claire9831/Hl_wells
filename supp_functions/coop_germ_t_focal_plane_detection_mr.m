% this code assumes a filled-out .csv file with optimal z-stack data exists

folderpath = 'Z:\2016-11-16';


directory = dir(folderpath);
directory = directory(3:end);


numericaldates = zeros(29,1);
for i_t = 1:29,
    directoryname = directory(i_t).name;
    wherearehyphens = regexp(directoryname,'-');
    directoryname(wherearehyphens)=[];
    numericaldates(i_t) = str2double(directoryname(1:12));
end;

[~,inds] = sort(numericaldates);
directory = directory(inds);

working_directory = dir([folderpath '\' directory(1).name]); % find directory for this time point
working_directory = working_directory(3:end);

filename = 'Z:\2016-11-16\2016-11-17-09-09-9717\singlewells\allwellsdata.csv';
wellxyfile = 'Z:\2016-11-16\2016-11-17-09-09-9717\singlewells\allwellsxy.mat';
load(wellxyfile); % this is where the data x_all and y_all come from

% ---------------------- Input data here -----------------------------
well_begin = 101; % beginning of range of wells to be analyzed
well_end = 200; % end of range of wells to be analyzed
% --------------------------------------------------------------------

optimal_z = csvread(filename,well_begin,1,['B' num2str(well_begin+1, '%i') '..B' num2str(well_end+1,'%i')]); % locate all of the data from this
num_wells = size(optimal_z,1); % the number of wells that have been manually analyzed
trackedwellinds = find(optimal_z>0); % the list of all wells that need to be tracked

% changed 28 to 29, the number of time points
best_z = zeros(num_wells,29); % will store the optimal focal plane at each time point

% maybe throw this out; I redefined the size here to deal with further
% calls to optimal_z
% optimal_z = optimal_z(1:num_wells);
best_z(:,1) = optimal_z; % initialize with the first time point

max_z_drift = 2; %  the max number of focal planes the optimal focus can shift by between frames



% first go through every single image and extract all of the wells that are needed
tic;

for i_t = 2:29,
    mkdir([folderpath '\' directory(i_t).name '\singlewells']);
    good_z = best_z(:,i_t-1);
    good_z = good_z(good_z > 0); % removing bad well z positions from consideration
    whichz = (min(good_z)-max_z_drift):(max(good_z)+max_z_drift);
    whichz = whichz(whichz > 0); % removes invalid z positions from the range
    % only analyze z-planes if they are within the dynamic range of at
    % least one well
    working_directory = dir([folderpath '\' directory(i_t).name]); % find directory for this time point
    working_directory = working_directory(3:end);
    
    for i_z = whichz
        % i_t+1 not needed; changed loop from 2:28 to 2:29
        
        % load image containing all panels at this z-plane
        im_thisz = imread([folderpath '\' directory(i_t).name '\' directory(i_t).name '_z' num2str(i_z,'%02i') '_ORG.tif']); % read in the image containing all panels, for this z and t value
        
        is_this_z = find(abs(best_z(:,i_t-1)-i_z)<=max_z_drift); % which wells could drift into this z-plane?
        % the below variable appears to be unnecessary
        % whichwells_inds = find(is_this_z); % indices of all tracked wells that could drift into this z-plane
        whichpanel_x = zeros(num_wells,1);
        whichpanel_y = zeros(num_wells,1);
        
        whichpanel_x(is_this_z) = floor(x_all(is_this_z+well_begin-1)/1388) + 1; % locate which panel each tracked well is in
        whichpanel_y(is_this_z) = floor(y_all(is_this_z+well_begin-1)/1040) + 1;
        
        whichpanel_ind = whichpanel_y+10*(whichpanel_x-1); % linear indices of panels containing wells
        whichpanel_ind(whichpanel_x==0) = 0;
        
        allpanelstoscan = unique(whichpanel_ind); % all panels that need to be analyzed
        allpanelstoscan(allpanelstoscan==0) = [];
        
        for panel_ind = allpanelstoscan.', % loop through panels (based on linear index)
            [panel_y,panel_x]=ind2sub([10 10],panel_ind); % find x-y locations of each panel
            im_thispanel = im_thisz((panel_y-1)*1040+1:panel_y*1040, ...
                (panel_x-1)*1388+1:panel_x*1388); % isolate sub-image for this panel
            
            % find wells that are in this focal range and in this panel
            
            well_centers = false(1040,1388);
            
            [x_thispanel,y_thispanel] = coopgerm_wellregister(im_thispanel);
            well_centers(sub2ind(size(well_centers),ceil(y_thispanel),ceil(x_thispanel))) = true; % has 1 only at the detected centers of wells
            wellsthispanel_bw = imdilate(well_centers,strel('square',220));
            
            % changed wells_thispanelbw to wellsthispanel_bw
            wellsthispanel_lab = bwlabel(wellsthispanel_bw);
            
            % changed i_panel to panel_ind
            is_wellinthispanel = whichpanel_ind==panel_ind; % find the wells (from the reduced, in-focus set) that are contained in this panel
            
            
            for i_well = find(is_wellinthispanel)', % changed call of x_all, y_all to be the i_wellth entry, should have the desired result.
                % i_well will always loop over the desired subset of
                % is_this_z, and so can always be called by x_all
                x_local = x_all(i_well+well_begin-1)-(panel_x-1)*1388; % find where well is located in the panel (using local coordinate system)
                y_local = y_all(i_well+well_begin-1)-(panel_y-1)*1040;
                
                well_lab = wellsthispanel_lab(y_local,x_local); % find which well in this panel matches the reference well
                wellsthislab = (wellsthispanel_lab == well_lab); % find all pixels that belong to this labeled well
                
                props = regionprops(wellsthislab,'Centroid');
                xy_thiswell = [props.Centroid];
                
                x_local = xy_thiswell(1);  y_local = xy_thiswell(2);
                
                x_all(i_well+well_begin-1) = ceil(x_local+(panel_x-1)*1388);
                y_all(i_well+well_begin-1) = ceil(y_local+(panel_y-1)*1040);
                
                isrowinthiswell = any(wellsthislab,2); % find y and x extent of the matching well
                Nrows = sum(isrowinthiswell);
                
                
                iscolinthiswell = any(wellsthislab,1);
                Ncols = sum(iscolinthiswell);
                
                % creates a padded well image if a well is detected,
                % otherwise it creates a properly-sized square of zeros
                if well_lab ~= 0
                    imthiswell = im_thispanel(isrowinthiswell,iscolinthiswell);
                    
                    
                    if Nrows<220 && isrowinthiswell(1),
                        imthiswell = [zeros(220-Nrows,size(imthiswell,2),'uint16'); imthiswell];
                    elseif Nrows<220 && isrowinthiswell(1040),
                        imthiswell = [imthiswell; zeros(220-Nrows,size(imthiswell,2),'uint16')];
                    end;
                    
                    if Ncols<220 && iscolinthiswell(1),
                        imthiswell = [zeros(220,220-Ncols,'uint16') imthiswell];
                    elseif Ncols<220 && iscolinthiswell(1388),
                        imthiswell = [imthiswell zeros(220,220-Ncols,'uint16')];
                    end
                    
                else
                    imthiswell = zeros(220,220,'uint16');
                end
                
                
                
                
                
                
                imwellfile = [folderpath '\' directory(i_t).name '\singlewells\well' num2str(i_well+well_begin-1,'%04i') '.tif'];
                if i_z==max(1,best_z(i_well,i_t-1)-max_z_drift),
                    imwrite(imthiswell,imwellfile,'WriteMode','Overwrite');
                else
                    imwrite(imthiswell,imwellfile,'WriteMode','Append');
                end;
                
            end;
            
        end;
        
    end; % now we have built image stacks covering the dynamical range of z for all tracked wells at this specific i_T
    
    % once we have looped over z planes, loop over all wells that are to be
    % tracked, and try to find which single image matches best to the
    % current reference image
    
    for i_well = trackedwellinds',
        
        %         try
        if i_t==2,
            imwellreffile = [folderpath '\' directory(1).name '\singlewells\well' num2str(i_well+well_begin-1,'%04i') '.tif'];
            well_ref = double(imread(imwellreffile,best_z(i_well,1)));
            
        else
            imwellreffile = [imwellfocuseddir 'well' num2str(i_well+well_begin-1,'%04i') '.tif'];
            well_ref = double(imread(imwellreffile,i_t-1));
        end;
        
        well_ref(well_ref==0) = NaN;
        %         catch
        %             keyboard;
        %         end;
        % go through each image from the stack of candidate in focus images
        
        %         origcorrcoeff = corrcoeff;
        well_refav = nanmean(well_ref(:));
        well_refstd = nanstd(well_ref(:));
        
        %         well_origrefav = mean(well_origref(:));
        %         well_origrefstd = std(well_origref(:));
        
        i_z = 1:2*max_z_drift+1;
        z=best_z(i_well,i_t-1)+i_z-max_z_drift-1;
        z(z<1)=1;
        z(z>31)=31;
        z = unique(z);
        corrcoeff = zeros(length(z),1);
        for i_z = 1:length(z),
            imwellfile = [folderpath '\' directory(i_t).name '\singlewells\well' num2str(i_well+well_begin-1,'%04i') '.tif'];
            
            imwell = double(imread(imwellfile,i_z));
            
            % if the well is all 0, then corrcoeff is set to a very large negative number, 
            % otherwise the corrcoeff is computed
            if sum(imwell(:)) ~= 0
                imwell(imwell==0) = NaN;
                
                imwell_av = nanmean(imwell(:));
                imwell_std = nanstd(imwell(:));
                %             corrcoeff(i_z) = sum(sum((imwell-imwell_av))); % /imwell_std/well_refstd;
                corrcoeff(i_z) = nansum(nansum((imwell-imwell_av).*(well_ref-well_refav)))/imwell_std/well_refstd;
                %             origcorrcoeff(i_z) = sum(sum((imwell-imwell_av).*(well_origref-well_origrefav))); % /imwell_std/well_origrefstd;
            else 
                corrcoeff(i_z) = -1000;
            end
        end;
        [best_corr,best_corr_z] = max(corrcoeff);
        %         [~,best_corr_z] = max(origcorrcoeff);
        
        % if no correlation computation is performed, the well is not
        % written and best_z remains 0
        if best_corr ~= -1000
            best_z(i_well,i_t) = best_z(i_well,i_t-1)+best_corr_z-max_z_drift-1;
            if best_z(i_well,i_t) < 0
                best_z(i_well,i_t) = 0;
            end
            imwellfocuseddir = [folderpath '\singlewells\'];
            if ~exist(imwellfocuseddir,'dir'),
                mkdir(imwellfocuseddir);
            end;
            
            imwellfocusedfile = [imwellfocuseddir 'well' num2str(i_well+well_begin-1,'%04i') '.tif'];
            %         imwellorigfocusedfile = [imwellfocuseddir 'origzfocused-well' num2str(i_well+well_begin-1,'%04i') '.tif'];
            
            imwellfile = [folderpath '\' directory(i_t).name '\singlewells\well' num2str(i_well+well_begin-1,'%04i') '.tif'];
            imwellinfocus = imread(imwellfile,best_corr_z);
            %         imwellinorigfocus = imread(imwellfile,best_corr_z);
            
            if i_t==2,
                well_ref = uint16(well_ref);
                imwrite(well_ref,imwellfocusedfile,'WriteMode','OverWrite');
                %             imwrite(well_ref,imwellorigfocusedfile,'WriteMode','OverWrite');
                
            end;
            % should these lines be in an else statement?
            imwrite(imwellinfocus,imwellfocusedfile,'WriteMode','Append');
            %         imwrite(imwellinorigfocus,imwellorigfocusedfile,'WriteMode','Append');
        end
        % else best_z in this position is 0 by initialization of best_z, and is then subsequently ignored later
    end
    
    
end;
t_alg = toc;
t_alg = t_alg / 3600;

