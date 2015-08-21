function [flag, timeStamp,newFilename,initializeFolder] = zenPoller_5Frame(PollDir,initializeFolder,expectedFileBase,counter)


%counter = counter_I0
if counter == 1
    expectedFileName = [expectedFileBase '.lsm'];
else
    expectedFileName = [expectedFileBase, num2str(counter-2) '.lsm'];
end

listing = dir([PollDir, '\' expectedFileName ]); % load the current folder
formatOut = 'HH:MM:SS'; % set format for timestamp can be usefull for consitency checks

if isempty(listing) ==1
                flag = 0;
                lastElement = [];
                timeStamp = [];
                newFilename = [];
                initializeFolder = dir(PollDir);

elseif listing.name == expectedFileName

                flag = 1;
                newFilename =  listing.name;

                timeStamp = datestr(listing.datenum,formatOut);
                initializeFolder = dir(PollDir);
else
                flag = 0;
                lastElement = [];
                timeStamp = [];
                newFilename = [];
                initializeFolder = dir(PollDir);

end


end







  
 
  