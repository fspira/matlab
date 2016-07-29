%%%%%%% This program aligns linescans - the workspace need to be loaded


 [mm pp nn] = size(greenStack);
 [greenOut redOut midOut] = doAverageFlanks(redMax1,redMax2,yGreen1,yGreen2,yRed1,yRed2)


 %%%%%% detect maxium size
 for lauf=1:nn
  
     lengthStore(lauf) = length(greenOut{lauf})
 
 end
 
 maxVector = max(lengthStore)
 
 %%%%%% sort linescans
 
 centerPoint = max(midOut);

 midOut
 
 centerMatrix = zeros(nn,maxVector)
 clear greenStore
 for lauf =1:nn
    
     
     %%%% fill zeros to the left side
     centerOffset = centerPoint - midOut(lauf)
     zerosFill = zeros(1,centerOffset);
     
     greenOutTmp = greenOut{lauf};
     redOutTmp = redOut{lauf};
     
     greenOutCat = cat(1,zerosFill',greenOutTmp)
     redOutCat = cat(1,zerosFill', redOutTmp)
    
     %%%% fill zeros to the right side
     centerOffset = maxVector - length(greenOutCat)
     if centerOffset < 0
         greenOutCat = greenOutCat(1:maxVector)
         redOutCat = redOutCat(1:maxVector)
     else
     zerosFill = zeros(1,centerOffset);
     
     
     greenOutCat = cat(1,greenOutCat, zerosFill')
     redOutCat = cat(1,redOutCat, zerosFill')
     end
     
     greenStore(:,lauf) = greenOutCat
 
     redStore(:,lauf) = redOutCat
 
 end
 
 
 
 for lauf =1:nn
 
     plot(greenStore(:,lauf),'g')
     hold on
     plot(redStore(:,lauf),'r')
     pause(0.2)
 end
 
 close all
 
 
 figure(1)
 imshow(greenStore,[])
 colormap(jet)
 
 figure(2)
 
 imshow(redStore,[])
 colormap(jet)
 
 %open greenStore
 relDistance = ((1:maxVector)*voxelX_mum) - (maxVector*voxelX_mum)/2
relDistance = relDistance'
open relDistance
timeVec = timeVec';

 