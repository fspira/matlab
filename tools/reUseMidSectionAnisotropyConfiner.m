FurrowRawMerge = (Furrow1Raw_S1_Store +Furrow1Raw_S2_Store)/2
PoleRawMerge = (Pole1Raw_S1_Store +Pole1Raw_S2_Store)/2

plot(FurrowRawMerge,'b')
hold on
plot(PoleRawMerge,'r')
axis([0 30 0 65000])
%StressVer = 1.312822185;
%StressHor = 0.541647473;

%StressAnisoVer = 0.092815147;
%StressAnisoHor = -0.178383408;


%PlasticCorrFactor = 1.32


%%%%%%% Calib M21VentralIngression 160421
StressAnisoHor =  0.428556937	
StressAnisoVer = 1.220190516
PlasticCorrFactor = 1.321237169



for lauf = 1:p
    %%%%% Calculate anisotropy
    
    
     anisoFurrow1(lauf) = (Furrow1Raw_S1_Store(lauf) - (PlasticCorrFactor * Furrow1Raw_S2_Store(lauf))) ...
         ./ (Furrow1Raw_S1_Store(lauf) + (2*PlasticCorrFactor*Furrow1Raw_S2_Store(lauf)))
     
     
     anisoFurrow2(lauf) = (Furrow2Raw_S1_Store(lauf) - (PlasticCorrFactor * Furrow2Raw_S2_Store(lauf))) ...
         ./ (Furrow2Raw_S1_Store(lauf) + (2*PlasticCorrFactor*Furrow2Raw_S2_Store(lauf)))
   
    anisoPole1(lauf) = (Pole1Raw_S1_Store(lauf) - (PlasticCorrFactor * Pole1Raw_S2_Store(lauf))) ...
         ./ (Pole1Raw_S1_Store(lauf) + (2*PlasticCorrFactor*Pole1Raw_S2_Store(lauf)))
     
    
    anisoPole2(lauf) =  (Pole2Raw_S1_Store(lauf) - (PlasticCorrFactor * Pole2Raw_S2_Store(lauf))) ...
         ./ (Pole2Raw_S1_Store(lauf) + (2*PlasticCorrFactor*Pole2Raw_S2_Store(lauf)))
     
    

     
end

%%%%% Normalize anisotropy
anisoFurrowMerge = (anisoFurrow1+ anisoFurrow2)/2;
anisoPoleMerge = (anisoPole1+anisoPole2)/2;

anisoFurrowNorm = (anisoFurrowMerge - StressAnisoHor) ./ (StressAnisoVer - StressAnisoHor)
anisoPoleNorm = (anisoPoleMerge - StressAnisoHor) ./ (StressAnisoVer - StressAnisoHor)

      saveVariables = {};

            saveVariables{1} = timeVec';%time{1};
            saveVariables{2} = contractileRingDistanceTmp';
         
            saveVariables{3} = FurrowRawMerge';
            saveVariables{4} = PoleRawMerge';
            
            
            
            saveVariables{5} = anisoFurrowMerge';%time{1};
            saveVariables{6} = anisoPoleMerge';
         
            saveVariables{7} = anisoFurrowNorm';
            saveVariables{8} = anisoPoleNorm';
     
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, contractile_ring_distance, RAWFurrow, RAWPole,',...
                'anisoFurrowRAW, anisoPoleRAW, anisoFurrowNorm, anisoPoleNorm'];
            
            outid = fopen([tifFilename,'Box_RatioAnalysis_midRAW_II.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            
            dlmwrite ([tifFilename,'Box_RatioAnalysis_midRAW.csv'],csvData,'roffset',1,'-append')
            
            
              save([tifFilename,'_aniso.mat'])
            