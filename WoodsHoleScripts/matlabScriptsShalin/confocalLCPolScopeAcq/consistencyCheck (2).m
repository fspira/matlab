function [consistencyFlag] = consistencyCheck(PollDir, expectedFileBase, counter)
  % compare "old" and "new" folder, this function will test for changes in
  % the folder structure during LC-change
  %counter = counter_I0
  if counter == 1
      consistencyFlag = 1;
  elseif counter == 2
     
      expectedFileName = [expectedFileBase, num2str(counter-1) '.lsm'];
      listingOld = dir([PollDir, '\' expectedFileName ]); % load the current folder
      if isempty(listingOld) == 1
          consistencyFlag = 1;
          
      else
          consistencyFlag = 0;
          
      end
      
  elseif counter > 2
     
      expectedFileName = [expectedFileBase, num2str(counter-1) '.lsm'];
      
      listingOld = dir([PollDir, '\'  expectedFileName ]); % load the current folder
      if isempty(listingOld) == 1
          consistencyFlag = 1;
          
      else
          consistencyFlag = 0;
          
      end
  end
  
  
  
end
