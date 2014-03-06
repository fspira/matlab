%%%%% Normalize ratio to metaphase  


   

        FurrowRawMeanStore_S1 = (Furrow1Raw_S1_Store + Furrow2Raw_S1_Store)./2
        FurrowRawMeanStore_S2 = (Furrow1Raw_S2_Store +  Furrow2Raw_S2_Store)./2
        
       
     
    
       meanPoleS1Raw =  (Pole1Raw_S1_Store +  Pole2Raw_S1_Store)./2
       meanPoleS2Raw =  (Pole1Raw_S2_Store +  Pole2Raw_S2_Store)./2
       
       timeZero = find(timeVec == 0)
       
       
       FurrowS1Norm = FurrowRawMeanStore_S1 ./  meanPoleS1Raw(timeZero);
       FurrowS2Norm = FurrowRawMeanStore_S2 ./ meanPoleS2Raw(timeZero);
       
       
       PoleS1Norm =  meanPoleS1Raw ./ meanPoleS1Raw(timeZero);
       PoleS2Norm =  meanPoleS2Raw ./ meanPoleS2Raw(timeZero);
        
       ratioFurrowNorm = FurrowS1Norm ./ FurrowS2Norm;
       
       ratioPoleNorm =   meanPoleS1Raw ./ meanPoleS2Raw 
      
        
       plot(timeVec,FurrowS1Norm,'r')
       hold on
       plot(timeVec, PoleS2Norm,'b')

       
       plot(timeVec, FurrowS1Norm,'r')
       hold on
       plot(timeVec, PoleS1Norm,'b')

