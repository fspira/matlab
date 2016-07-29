
function [expectedValues,flank1DiffFromExpected,alphaDegNew] =  doCalculateExpectedIntenisties_splineFit5Knots_Div_M21(funCoords,ratioNormFlank1Tmp,fitresult,testImg,maxMean)

lauf = 1;
   
    
    uNew = [];
    vNew = [];
    alphaDeg1New = [];
    correctedFlank1 = [];
    flank1DiffFromExpected = [];
     expectedValues = [];
    
    
  % funCoords =  flank2Tmp
   
 varLength =  length(funCoords);


%%%%%%%%%
%%%%%%%%% Determine zeroLength - this is important to determine the number
%%%%%%%%% of knots
%%%%%%%%%
   zeroLengthPara =0;
   subtractor = 0;
  while zeroLengthPara == 0
      
      paraTest = (varLength - subtractor) / 5;
      paraTestEval = round(paraTest/5)*5;
      
      if paraTest == paraTestEval
          zeroLengthPara = 1;
      end
      
      subtractor = subtractor+1;

  end
  
  
  zeroLength = varLength - subtractor+1;
    %zeroLength = 540;

    cv2_StoreSub(1:2,1:zeroLength)  = zeros;

    cdv2_StoreSub(1:2,1:zeroLength)  = zeros;
    skelfitKnotsStore2Sub(1,1:zeroLength) = zeros;
    
 
marker = 1;
 subrun = 1;
 for loopRun = 1:5
    
     testImg(:,:) = 0;

     for lauf =marker:zeroLength-5+marker

         testImg(funCoords(1,lauf),funCoords(2,lauf)) = 255;

     end

 BWTest = im2bw(testImg);
 I3 = bwmorph(BWTest,'bridge',8);

 %imshow(I3,[])
 

     %%%% Determine the number of knots - interpolate over 5 pixels

     knotSize = round(zeroLength/5);

     [skelfit, skeldist, rowI, colI, rowseed, colseed, skelI]=splinefitskel(I3,knotSize,5);

    fnplt(skelfit,'or',1); hold on
    cv = fnval(skelfit, skelfit.knots);
    cdv = fnval(fnder(skelfit), skelfit.knots);
    quiver(cv(1,:),cv(2,:), cdv(1,:),cdv(2,:));


    cv2_StoreSub(1:2,marker:5:end) = cv(1:2,1:knotSize);
    cdv2_StoreSub(1:2,marker:5:end) = cdv(1:2,1:knotSize);
    skelfitKnotsStore2Sub(1,marker:5:end) = skelfit.knots(1:knotSize);
          

    subrun = subrun+2;
    marker = marker+1;

 end
    quiver(cv2_StoreSub(1,:),cv2_StoreSub(2,:), cdv2_StoreSub(1,:),cdv2_StoreSub(2,:));hold on
    
    
    %%%%%%%%% Esay angle detection method by simply calcualting the angle
    %%%%%%%%% between consecuteive points
    
  for subrun= 1:length(cdv2_StoreSub(1,:))
      
     
      
  % v = sqrt((cdv2_StoreSub(1,subrun))^2 + cdv2_StoreSub(2,subrun)^2);

   % alphaDegNew(subrun) = acosd(cdv2_StoreSub(1,subrun)/v);
   % IntenistyFlank1(subrun) = greenOverlayStoreNorm(round(cv(1,subrun)),round(cv(2,subrun)))
    rTmp = mod(atan2(cdv2_StoreSub(2,subrun),cdv2_StoreSub(1,subrun)),2.*pi);
    %if rTmp < 0
   % rTmp = atan2(cdv2_StoreSub(2,subrun),cdv2_StoreSub(1,subrun))/pi*180 + 180;
   % end
   rOut( subrun)= rTmp;
    
   
    
   alphaDegNew(subrun) = rTmp*180/pi;
     
  end
  
 
anisoMetaCorr = feval(fitresult,30);

  splineEnd = length(ratioNormFlank1Tmp);
  
  if zeroLength > length(ratioNormFlank1Tmp)
      
     flankCrop1 = ratioNormFlank1Tmp(1:length(ratioNormFlank1Tmp));
     zeroLength = length(flankCrop1);
  else
      
      flankCrop1 = ratioNormFlank1Tmp(1:zeroLength);
  end
      
  
  if length(cdv2_StoreSub) < zeroLength
      
      zeroLength = length(cdv2_StoreSub);
      
  end
    
  
  
        for subrun = 1:zeroLength%length(flank1)


        expectedValues(subrun) =  maxMean - (feval(fitresult,(alphaDegNew(subrun))));
      %  flank1DiffFromExpected(subrun) =  flankCrop1(subrun) + (maxMean - (feval(fitresult, alphaDegNew(subrun))));
         
        
         
      anisoCorrectionFactor = (feval(fitresult, (alphaDegNew(subrun))-10));
    
   if  flankCrop1(subrun) >= (feval(fitresult, (alphaDegNew(subrun))-10))
  
 
    
        flank1DiffFromExpected(subrun) =  flankCrop1(subrun) ./ anisoCorrectionFactor;
    else
       flank1DiffFromExpected(subrun) =  flankCrop1(subrun) ./ anisoCorrectionFactor;
    end
        
        
      %   flank1DiffFromExpectedDivStress(subrun) = 1-( (flankCrop1(subrun)+ (maxMean - (feval(fitresult, alphaDegNew(subrun)))))./ (feval(fitresultstressFibers,alphaDegNew(subrun))));
       
       % flank1DiffFromExpectedDiv(subrun) =  flankCrop1(subrun) ./(feval(fitresult, alphaDegNew(subrun)));
   
        end
        
        
end
        
     
        

