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

filename = [folderpath '\' directory(1).name '\singlewells\allwellsdata.csv'];
wellxyfile = [folderpath '\' directory(1).name '\singlewells\allwellsxy.mat'];
successfile = [folderpath '\' directory(1).name '\singlewells\success_flag.mat'];
load(wellxyfile); % this is where the data x_all and y_all come from
load(successfile);

wellstobeanalyzed = success_flag < 20; % MAKE SURE TO UPDATE THIS AT END OF CODE
well_begin = 1; 


optimal_z = csvread(filename,1,1); % locate all of the data from this
optimal_z = optimal_z(:,1);
num_wells = size(optimal_z,1); % the number of wells that have been manually analyzed
optimal_z(~wellstobeanalyzed) = 0;
trackedwellinds = find(optimal_z>0); % the list of all wells that need to be tracked

% changed 28 to 29, the number of time points
best_z = zeros(num_wells,29); % will store the optimal focal plane at each time point

% optimal_z = optimal_z(1:num_wells);
best_z(:,1) = optimal_z; % initialize with the first time point