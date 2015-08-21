function [consistencyFlag] = consistencyCheck(PollDir,initializeFolder)
  % compare "old" and "new" folder, this function will test for changes in
  % the folder structure during LC-change
 
if isequal(initializeFolder,Polldir) == 1
    consistencyFlag = 1;
else
    consistencyFlag = 0;
    
end
