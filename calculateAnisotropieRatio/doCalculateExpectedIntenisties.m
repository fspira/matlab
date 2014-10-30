function [expectedValues,flank1DiffFromExpected,alphaDeg1New] =  doCalculateExpectedIntenisties(I1_coords,ratioNormFlank1Tmp,fitresult,flankCrop1)

lauf = 1
    
    
    uNew = []
    vNew = []
    alphaDeg1New = [];
    correctedFlank1 = []
    flank1DiffFromExpected = [];
     expectedValues = [];
    
    
    
    
   
   
    
    %%%%%%%%% Esay angle detection method by simply calcualting the angle
    %%%%%%%%% between consecuteive points
    
    for subrun = 1:(length(I1_coords))-5
        
        uNew(subrun) = I1_coords(1,subrun+5) - I1_coords(1,subrun);
        vNew(subrun) = I1_coords(2,subrun+5) - I1_coords(2,subrun);
    
    end
    
    
    
   % quiver(I1_coords(1,1:end-5),I1_coords(2,1:end-5),  uNew(1,:),vNew(1,:),2);hold on
    
    %spline = cscvn(I1_coords);
    %fnplt(spline,'or',2.5); hold on

    
    %%%%%% Calculate the angle between the individual knots
    
    for subrun = 1:length(vNew)
    
            v1New = sqrt((vNew(1,subrun))^2 + uNew(1,subrun)^2);
            alphaDeg1New(subrun) = acosd(uNew(1,subrun)/v1New)
                
    end
    
    
    
    
      
      flankCrop1 = ratioNormFlank1Tmp(1:length(vNew(1,:)));
      
      
    
        for subrun = 1:length(vNew(1,:))%length(flank1)


        expectedValues(subrun) =  feval(fitresult,alphaDeg1New(subrun))
        flank1DiffFromExpected(subrun) =   (feval(fitresult, alphaDeg1New(subrun))) -flankCrop1(subrun)

subrun
        end
        
        
end
        
     
        

