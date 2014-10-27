function  [greenOverlayStore]= doGenerateAnisotropieOverlay(greenNorm, ratioFurrow1, ratioFurrow2,flank1Store,flank2Store,AnalysisEnd)

 
 
 for subrun = 1:AnalysisEnd
 
     greenOverlay = greenNorm(:,:,1);
     greenOverlay(:,:) = 0;
     greenOverlay = double(greenOverlay);


     lauf = 1


     ratioFurrow1Tmp = ratioFurrow1{subrun}
     ratioFurrow2Tmp = ratioFurrow2{subrun}
     flank1Tmp = flank1Store{subrun};
     flank2Tmp = flank2Store{subrun};


     for lauf = 1: length(flank1Tmp)-1

         greenOverlay(flank1Tmp(lauf,2),flank1Tmp(lauf,1)) = ratioFurrow1Tmp(lauf);
         greenOverlay(flank2Tmp(lauf,2),flank2Tmp(lauf,1)) = ratioFurrow2Tmp(lauf);
         
         %pause(0.2)
        greenOverlayStore(:,:,subrun) = greenOverlay; 

     end
    
 end
 
 
 
 
 imshow(greenOverlayStore(:,:,11),[]);
         colormap(jet)
       