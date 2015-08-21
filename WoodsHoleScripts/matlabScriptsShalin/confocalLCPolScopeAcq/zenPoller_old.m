function [flag, timeStamp,newFilename,initializeFolder] = zenPoller(PollDir,initializeFolder)

listing = dir(PollDir); % load the current folder
formatOut = 'HH:MM:SS'; % set format for timestamp can be usefull for consitency checks

% test for changes in the folder structure, by comparing the "old"
% initialized folder with the current "new" state of the folder

if isequal(listing,initializeFolder) 
    flag = 0; % if old and new folder are the same, do nothing
    lastElement = [];
    timeStamp = [];
    newFilename = [];
    initializeFolder = dir(PollDir);
else
    
    % get the name of the last element in the folder, since zen
    % nomenclature is incremental, the last element in the folder will be
    % the latest file.
    
    lastElementOld = size(initializeFolder); 
    lastElement = size(listing);
    oldFilename = initializeFolder(lastElementOld(1)).name;
    newFilename = listing(lastElement(1)).name;
    
    
    
    %disp(newFilename)
    % check for file length, if file length is > 2 the folder will
    % be empty
    if  length(newFilename) > 2  
        %disp('flag1')
        % compare the names of last files in the "old" and "new" folder if
        % both are the same do nothing this is to test whether the change
        % was indeed a real incremental ZEN file, or some temporary file
        % (which would not appear in the last position of the folder list)
        
        if strcmp(newFilename,oldFilename)
            %disp('flag2')
            flag = 0;
            lastElement = [];
            timeStamp = [];
            newFilename = [];
            initializeFolder = dir(PollDir);
            % if the new file is different from the older file, trigger the
            % change flag and set the "old" folder to "new" folder
            
        else
            %disp('flag3')
            flag = 1;
            newFilename = listing(lastElement(1)).name;
            timeStamp = datestr(listing(lastElement(1)).datenum,formatOut);
            initializeFolder = dir(PollDir);
        end
    else
        % disp('flag4')
        % if the floder is empty, do nothing
        flag = 0;
        lastElement = [];
        timeStamp = [];
        newFilename = [];
        initializeFolder = dir(PollDir);
        
    end
    
end




  
  
  
 
  