% M. Roper version of image making file. Loops through images and isolates
% single wells, stacking images taken at different z-positions, so we can
% identify the optimal focus for each well
%
% introduce 

folderpath = 'Z:\2017-01-18';
grid_size = [10 10]; % number of x-y- panels that are combined to make the total image
numRows = 1040; % size of individual panels
numCols = 1388;

alldirectories = dir(fullfile(folderpath, '2017*'));
N_t = numel(alldirectories);
allc1_files = dir(fullfile([folderpath '\' alldirectories(1).name],'*c1_ORG.tif'));
if numel(allc1_files)==0,
    all_files = dir(fullfile([folderpath '\' alldirectories(1).name],'*_ORG.tif'));
    N_z = numel(all_files);
    is_2channel = false;
else
    N_z = numel(allc1_files);
    is_2channel = true;
end;

stack = cell(N_z,N_t);
well_squares = cell(N_z,N_t);

directory = dir(folderpath);
directory = directory(3:end);


numericaldates = zeros(N_t,1);
for i_t = 1:N_t,
    directoryname = directory(i_t).name;
    wherearehyphens = regexp(directoryname,'-');
    directoryname(wherearehyphens)=[];
    numericaldates(i_t) = str2double(directoryname(1:12));
end;

[~,inds] = sort(numericaldates);
directory = directory(inds);

% new structure: 1. Find all wells for the first time point, and store
% images of them

for t = 1, 
    singlewelldir = [folderpath '\' directory(t).name '\singlewells'];
    if ~exist(singlewelldir,'dir'),
        mkdir(singlewelldir);
    end;
    
    if is_2channel,
        working_directory = dir(fullfile([folderpath '\' directory(t).name],'*c2_ORG.tif')); % find directory for this time point
        gfpfiles = dir(fullfile([folderpath '\' directory(t).name],'*c1_ORG.tif')); % find directory for this time point
    else
        working_directory = dir([folderpath '\' directory(t).name]); % find directory for this time point
        working_directory = working_directory(3:end);
    end;
        N_wells_z0 = zeros(10,10); % number of wells that were detected in each panel in the first frame of the image
    x_all = zeros(100*16,1); y_all = zeros(100*16,1); % stores all the wells that were detected at z0
    
    for z = 1:N_z,
        im = imread([folderpath '\' directory(t).name '\' working_directory(z).name]); % find image corresponding to a single z-position
        if is_2channel,
            imgfp = imread([folderpath '\' directory(t).name '\' gfpfiles(z).name]); % find image corresponding to a single z-position
        end;
        N_wells = 0;
        
        for i_panel = 1:100,
            [i_y,i_x] = ind2sub([10 10],i_panel);
            im_panel = im((i_y-1)*numRows+(1:numRows),(i_x-1)*numCols+(1:numCols)); % we will find the wells in each panel individually, neglecting wells that span two panels
            if is_2channel,
                imgfp_panel= imgfp((i_y-1)*numRows+(1:numRows),(i_x-1)*numCols+(1:numCols)); % we will find the wells in each panel individually, neglecting wells that span two panels
            end;
            [x,y] = coopgerm_wellregister(im_panel);% find all of the wells within this frame
            
            x = round(x); y= round(y);
            
            squaresize = 220;
            marginsize = 0.5*squaresize;
            
            if z==1,
                
                is_ok = x>marginsize & x<=numCols-marginsize & y>marginsize & y<=numRows-marginsize; % only keep wells that are fully contained within this panel
                x = x(is_ok); y = y(is_ok); % clean up wells that have been detected
                
                N_wells_thispanel = length(x); % find the number of wells in this panel
                N_wells_z0(i_y,i_x) = N_wells_thispanel;
                x_all(N_wells+(1:N_wells_thispanel)) = (i_x-1)*numCols+x; % master list of all wells that were detected across all panels
                y_all(N_wells+(1:N_wells_thispanel)) = (i_y-1)*numRows+y;
                
            else % only worry about wells whose centers are outside of the image
                is_ok = x>=1 & x<=numCols & y>=1 & y<=numRows;
                x = x(is_ok); y = y(is_ok); % clean up wells that have been detected
            end;
            
            inds_ok = sub2ind([numRows numCols],y,x);
            
            wells_bw = false(size(im_panel));
            wells_bw(inds_ok) = 1;
            wells_bw = imdilate(wells_bw,strel('square',220));
            wells_label = bwlabel(wells_bw);
            
            if z==1,
                for i_well = 1:N_wells_thispanel;
                    imfilename = [folderpath '\' directory(t).name '\singlewells\well' num2str(N_wells+i_well,'%04i') '.tif'];
                    imvisfolder = [folderpath '\' directory(t).name '\singlewellsforvisualization'];
                    if ~exist(imvisfolder,'dir'),
                        mkdir(imvisfolder);
                    end;
                    imvisfilename = [imvisfolder '\well' num2str(N_wells+i_well,'%04i') '.tif'];
                    
                    thiswell_bw = (wells_label == wells_label(y(i_well),x(i_well))); % isolate the pixels in this well
                    isrowinthiswell = any(thiswell_bw,2);
                    
                    iscolinthiswell = any(thiswell_bw,1);
                    
                    Nrows = sum(isrowinthiswell);
                    Ncols = sum(iscolinthiswell);
                    
                    im_well = im_panel(isrowinthiswell,iscolinthiswell);
                    if is_2channel,
                        imgfp_well = imgfp_panel(isrowinthiswell,iscolinthiswell);
                    end;
                        
                    if Nrows<220 && isrowinthiswell(1),
                        
                        im_well = [zeros(220-Nrows,size(im_well,2),'uint16'); im_well];
                        if is_2channel,
                            imgfp_well = [zeros(220-Nrows,size(imgfp_well,2),'uint16'); imgfp_well];
                        end;
                        
                    elseif Nrows<220 && isrowinthiswell(1040),
                        
                        im_well = [im_well; zeros(220-Nrows,size(im_well,2),'uint16')];
                        if is_2channel,
                            imgfp_well = [imgfp_well; zeros(220-Nrows,size(imgfp_well,2),'uint16')];
                        end;
                    end;
                    
                    if Ncols<220 && iscolinthiswell(1),
                        
                        im_well = [zeros(220,220-Ncols,'uint16') im_well];
                        if is_2channel,
                            imgfp_well = [zeros(220,220-Ncols,'uint16') imgfp_well];
                        end;
                    elseif Ncols<220 && iscolinthiswell(1388),
                        
                        im_well = [im_well zeros(220,220-Ncols,'uint16')];
                        if is_2channel,
                            imgfp_well = [imgfp_well zeros(220,220-Ncols,'uint16')];
                        end;
                        
                    end;
                    
                    imwrite(im_well,imfilename);
                    im_well_adj = imadjust(im_well);
                    if is_2channel,
                        imgfp_well_adj = 50*imgfp_well;
                        im_well_adj = repmat([im_well_adj im_well_adj],[1 1 3]);
                        im_well_adj(1:220,221:end,2) = max(im_well_adj(1:220,221:end,2),imgfp_well_adj);
                    end;
                    
                    imwrite(im_well_adj,imvisfilename);
                                    end;
                N_wells = N_wells + N_wells_thispanel; % number of wells detected so far
            else % loop over all wells detected in z0, and find the matching well in z
                % add the image of the well in z to the image in z0
                for i_well = 1:N_wells_z0(i_y,i_x), % for each well in z0
                    x = x_all(N_wells+i_well)-(i_x-1)*numCols;
                    y = y_all(N_wells+i_well)-(i_y-1)*numRows;
                    
                    imfilename = [folderpath '\' directory(t).name '\singlewells\well' num2str(N_wells+i_well,'%04i') '.tif'];
                    imvisfilename = [folderpath '\' directory(t).name '\singlewellsforvisualization\well' num2str(N_wells+i_well,'%04i') '.tif'];

                    if wells_label(y,x)>0, % find which well matches this one in current z
                        thiswell_bw = (wells_label == wells_label(y,x)); % isolate the pixels in this well
                        isrowinthiswell = any(thiswell_bw,2);
                        iscolinthiswell = any(thiswell_bw,1);

                        if is_2channel,
                            imgfp_well = imgfp_panel(isrowinthiswell,iscolinthiswell);
                        end;
                        im_well = im_panel(isrowinthiswell,iscolinthiswell);
                    else % then for some reason, no well was detected in this frame
                        im_well = zeros(squaresize,'like',im_panel);
                        if is_2channel,
                            imgfp_well = zeros(squaresize,'like',im_panel);
                        end;
                        
                    end;

                    Nrows = sum(isrowinthiswell);
                    Ncols = sum(iscolinthiswell);
                        
                    if Nrows<220 && isrowinthiswell(1),
                        
                        im_well = [zeros(220-Nrows,size(im_well,2),'uint16'); im_well];
                        if is_2channel,
                            imgfp_well = [zeros(220-Nrows,size(imgfp_well,2),'uint16'); imgfp_well];
                        end;
                        
                        
                    elseif Nrows<220 && isrowinthiswell(1040),
                        
                        im_well = [im_well; zeros(220-Nrows,size(im_well,2),'uint16')];
                        if is_2channel,
                            imgfp_well = [imgfp_well; zeros(220-Nrows,size(imgfp_well,2),'uint16')];
                        end;
                    end;
                    
                    if Ncols<220 && iscolinthiswell(1),
                        
                        im_well = [zeros(220,220-Ncols,'uint16') im_well];
                        if is_2channel,
                            imgfp_well = [zeros(220,220-Ncols,'uint16') imgfp_well];
                        end;
                    elseif Ncols<220 && iscolinthiswell(1388),
                        
                        im_well = [im_well zeros(220,220-Ncols,'uint16')];
                        if is_2channel,
                            imgfp_well = [imgfp_well zeros(220,220-Ncols,'uint16')];
                        end;
                        
                    end;
                    
                    imwrite(im_well,imfilename,'Writemode','append');
                    im_well_adj = imadjust(im_well);
                    if is_2channel,
                        imgfp_well_adj = 50*imgfp_well;
                        im_well_adj = repmat([im_well_adj im_well_adj],[1 1 3]);
                        im_well_adj(1:220,221:end,2) = max(im_well_adj(1:220,221:end,2),imgfp_well_adj);
                    end;

                    imwrite(im_well_adj,imvisfilename,'WriteMode','append');
                    
                end;
                N_wells = N_wells + N_wells_z0(i_y,i_x); % number of wells that have been matched so far
            end;
            
        end;
        % once all panels have been looped over delete vacant spaces in the
        % arrays x_all and y_all
        if z==1,
            x_all(N_wells+1:end) = [];
            y_all(N_wells+1:end) = [];
            x_all = ceil(x_all);
            y_all = ceil(y_all);
            savefile = [folderpath '\' directory(t).name '\singlewells\allwellsxy.mat'];
            save(savefile,'x_all','y_all','N_wells_z0','N_t','N_z'); % store the coordinates of all the wells that were detected at z=1
            success_flag = zeros(length(x_all),1);
            successfile = [folderpath '\' directory(1).name '\singlewells\success_flag.mat'];
            save(successfile,'success_flag');
            
            
            xlsfile = [folderpath '\' directory(t).name '\singlewells\allwellsdata.csv'];
            if ~exist(xlsfile,'file'),
                fileID = fopen(xlsfile,'w');
                
                fprintf(fileID,'%s\n','well id,optimal z,N_spores1,area1,N_spores2,area2');
                
                for i_well=1:N_wells,
                    fprintf(fileID,'%s\n',['well' num2str(i_well,'%04i') ',,,,,']);
                end;
                
                fclose(fileID);
            end;
        end;
        
        
    end;
end;



