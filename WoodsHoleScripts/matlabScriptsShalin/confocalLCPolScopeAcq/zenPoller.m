function [flag, timeStamp,newFilename,initializeFolder,lastFileInFolder] = zenPoller(PollDir,initializeFolder,lastFileInFolder,fileBase)
subcounter = 1
subcounter = subcounter +1

listing = dir(PollDir) % load the current folder
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
    
    %lastElementOld = size(initializeFolder);
    %lastElement = size(listing);
    oldFilename = initializeFolder(end).name;
    newFilename = listing(end).name;
    
    
    
    
    %disp('flag1')
    % compare the names of last files in the "old" and "new" folder if
    % both are the same do nothing this is to test whether the change
    % was indeed a real incremental ZEN file, or some temporary file
    % (which would not appear in the last position of the folder list)
    
    %%% Identify number of the old file
    %%%
    
    %if oldFilename
    
    if strcmp(oldFilename,'.')
        
        fileCounterOld = -1
    elseif strcmp(oldFilename,'..')
        fileCounterOld = -1
    else
        suffixStrOld = findstr(oldFilename,'.lsm');
        prefixStrOld = findstr(oldFilename,fileBase(end));
        fileCounterOld = str2num(oldFilename(prefixStrOld+1:suffixStrOld-1));
        
        %%% Identify number of the new file
        
        suffixStrNew = findstr(newFilename,'.lsm');
        prefixStrNew = findstr(newFilename,fileBase(end));
        fileCounterNew = str2num(newFilename(prefixStrNew+1:suffixStrNew-1));
        
        %strcmp(newFilename,oldFilename)
        
        if strcmp(initializeFolder(end).name,listing(end).name)
            
            flag = 0; % if old and new folder are the same, do nothing
            lastElement = [];
            timeStamp = [];
            newFilename = [];
            initializeFolder = dir(PollDir);
            
            
            if  fileCounterNew > fileCounterOld
                
                % change flag and set the "old" folder to "new" folder
                %disp('flag3')
                flag = 1;
                newFilename = listing(end).name;
                timeStamp = datestr(listing(end).datenum,formatOut);
                initializeFolder = dir(PollDir);
                
                
            else
                
                %disp('flag2')
                flag = 0;
                lastElement = [];
                timeStamp = [];
                newFilename = [];
                initializeFolder = dir(PollDir);
                % if the new file is different from the older file, trigger the
                
            end
            
            
        end
    end
    
end

end






  
 
  