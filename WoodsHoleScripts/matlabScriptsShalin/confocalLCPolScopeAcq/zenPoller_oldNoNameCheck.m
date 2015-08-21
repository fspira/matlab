function [flag, timeStamp,newFilename,initializeFolder] = zenPoller_oldNoNameCheck(PollDir,initializeFolder)

listing = dir(PollDir); % load the current folder
formatOut = 'HH:MM:SS'; % set format for timestamp can be usefull for consitency checks

lastElement = size(listing);
newLastElement = listing(lastElement(1)).name;
% test for changes in the folder structure, by comparing the "old"
% initialized folder with the current "new" state of the folder

    if isequal(listing,initializeFolder)
        flag = 0; % if old and new folder are the same, do nothing
        lastElement = [];
        timeStamp = [];
        newFilename = [];
        initializeFolder = dir(PollDir);
    elseif length(newLastElement)<= 2

        flag = 0; % if old and new folder are the same, do nothing
        lastElement = [];
        timeStamp = [];
        newFilename = [];
        initializeFolder = dir(PollDir);

        % get the name of the last element in the folder, since zen
        % nomenclature is incremental, the last element in the folder will be
        % the latest file.


        %%% sort reference folder elements
    elseif length(listing) <= 8


        lastElementOld = size(initializeFolder);
        lastElement = size(listing);
        oldFilename = initializeFolder(lastElementOld(1)).name;
        newFilename = listing(lastElement(1)).name;

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

    elseif length(listing) > 8

            lengthDirectoryReference = length(initializeFolder)

            clear listDirectory  listDirectoryReference sortArrayReference
            clear  sortedListReference sortOutReference sortIdxReference

            for subrun = 1:length(initializeFolder)-3
                laufNameReference = initializeFolder(subrun+3).name
                listDirectoryReference{subrun} =  laufNameReference
                suffixStrSortReference = findstr( laufNameReference,'.lsm');
                prefixStrSortReference = findstr( laufNameReference,'T');


                sortArrayReference(subrun) = str2num(laufNameReference(prefixStrSortReference+1:suffixStrSortReference-1));
            end

            [sortOutReference sortIdxReference] = sort(sortArrayReference)

            sortedListReference{1} = initializeFolder(3).name

            for subrun = 1:length(sortIdxReference)
                sortedListReference{subrun+1} = listDirectoryReference{sortIdxReference(subrun)};
            end



            %%% sort monitored folder elements


            lengthDirectory = length(listing)

            clear listDirectory sortOut
            clear sortedList sortArray sortIdx

            for subrun = 1:length(listing)-3
                laufName = listing(subrun+3).name
                listDirectory{subrun} =  laufName
                suffixStrSort = findstr( laufName,'.lsm');
                prefixStrSort = findstr( laufName,'T');


                sortArray(subrun) = str2num(laufName(prefixStrSort+1:suffixStrSort-1));
            end

            [sortOut sortIdx] = sort(sortArray)

            sortedList{1} = listing(3).name

            for subrun = 1:length(sortIdx)
                sortedList{subrun+1} = listDirectory{sortIdx(subrun)};
            end



            oldFilename =  sortedListReference{end}
            newFilename =  sortedList{end}


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
                newFilename =  sortedList{end};

                timeStamp = datestr(listing(sortIdx(end)).datenum,formatOut);
                initializeFolder = dir(PollDir);
            end


    end
end







  
 
  