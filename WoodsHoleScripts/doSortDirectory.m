function [sortArray listDirectory] = doSortDirectory(fileList)
            for subrun = 1:length(fileList)-1
                laufName = fileList(subrun+1).name
                listDirectory{subrun} =  laufName
                suffixStrSort = findstr( laufName,'.lsm');
                prefixStrSort = findstr( laufName,'_');


                sortArray(subrun) = str2num(laufName(prefixStrSort+1:suffixStrSort-1));
            end
            
end