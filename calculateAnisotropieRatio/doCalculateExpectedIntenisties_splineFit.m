function [expectedValues,flank1DiffFromExpected,alphaDegNew] =  doCalculateExpectedIntenisties_splineFit(funCoords,ratioNormFlank1Tmp,fitresult,flankCrop1,testImg)

lauf = 1
   
    
    uNew = []
    vNew = []
    alphaDeg1New = [];
    correctedFlank1 = []
    flank1DiffFromExpected = [];
     expectedValues = [];
    
    
  % funCoords =  flank2Tmp
   
 varLength =  length(funCoords)


%%%%%%%%%
%%%%%%%%% Determine zeroLength - this is important to determine the number
%%%%%%%%% of knots
%%%%%%%%%
   zeroLengthPara =0;
   subtractor = 0;
  while zeroLengthPara == 0
      
      paraTest = (varLength - subtractor) / 6
      paraTestEval = round(paraTest/5)*5
      
      if paraTest == paraTestEval
          zeroLengthPara = 1
      end
      
      subtractor = subtractor+1

  end
  
  
  zeroLength = varLength - subtractor+1;
    %zeroLength = 540;

    cv2_StoreSub(1:2,1:zeroLength)  = zeros

    cdv2_StoreSub(1:2,1:zeroLength)  = zeros
    skelfitKnotsStore2Sub(1,1:zeroLength) = zeros
    
 
marker = 1
 subrun = 1
 for loopRun = 1:6
    
     testImg(:,:) = 0;

     for lauf =marker:zeroLength-5+marker

         testImg(funCoords(1,lauf),funCoords(2,lauf)) = 255;

     end

 BWTest = im2bw(testImg);
  I3 = bwmorph(BWTest,'bridge',8);

 imshow(I3,[])
 

     %%%% Determine the number of knots - interpolate over 5 pixels

     knotSize = round(zeroLength/6)

     [skelfit, skeldist, rowI, colI, rowseed, colseed, skelI]=splinefitskel(I3,knotSize,6)

    fnplt(skelfit,'or',1); hold on
    cv = fnval(skelfit, skelfit.knots);
    cdv = fnval(fnder(skelfit), skelfit.knots);
    quiver(cv(1,:),cv(2,:), cdv(1,:),cdv(2,:));


    cv2_StoreSub(1:2,marker:6:end) = cv(1:2,1:knotSize);
    cdv2_StoreSub(1:2,marker:6:end) = cdv(1:2,1:knotSize);
    skelfitKnotsStore2Sub(1,marker:6:end) = skelfit.knots(1:knotSize);
          

    subrun = subrun+2
    marker = marker+1

 end
   
    
    %%%%%%%%% Esay angle detection method by simply calcualting the angle
    %%%%%%%%% between consecuteive points
    
  for subrun= 1:length(cdv2_StoreSub(1,:))
      
     
      
   v = sqrt((cdv2_StoreSub(1,subrun))^2 + cdv2_StoreSub(2,subrun)^2);

    alphaDegNew(subrun) = acosd(cdv2_StoreSub(1,subrun)/v);
   % IntenistyFlank1(subrun) = greenOverlayStoreNorm(round(cv(1,subrun)),round(cv(2,subrun)))
     
  end

  splineEnd = length(ratioNormFlank1Tmp)
  
  if splineEnd > 300
      
     splineEnd = 300;
     
  else
      
      flankCrop1 = ratioNormFlank1Tmp(1:zeroLength);
  end
      
    
        for subrun = 1:zeroLength%length(flank1)


        expectedValues(subrun) =  feval(fitresult,alphaDegNew(subrun))
        flank1DiffFromExpected(subrun) =   (feval(fitresult, alphaDegNew(subrun))) -flankCrop1(subrun)

subrun
        end
        
        
end
        
     
        

