%%%% Plot the ROI used for analysis and save the file

tifFilename= [tifFilename,'_Box1_1'];

%anaOnset = 10;

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;


for lauf =1 :length(xFurrow1)

    greenStackRGB(:,:,:,lauf) = ind2rgb(S1(:,:,lauf),green);
    redStackRGB(:,:,:,lauf) = ind2rgb(S2(:,:,lauf),red);
             
    imgMergeTmp(:,:,:,lauf) = greenStackRGB(:,:,:,lauf) + redStackRGB(:,:,:,lauf);
    imgMerge(:,:,:,lauf) = normalizedImage3D(imgMerge(:,:,:,lauf));
    
    
end

  clear ratioFurrowStore ratioPoleStore  saveImg

for lauf = 1:length(xFurrow1)

    x = xFurrow1{lauf}
    y = yFurrow1{lauf}
    
   S1_In = imgMerge(:,:,:,lauf);   
   % S2_In = S2(:,:,:,lauf);   
    
 [BW_Box BWHor BWVer,imgMergeSaveVertical]=doDrawBoxAtCoordinates_Hor_1_1(x,y,S1_In,voxelX_mum);
 
 %%%%% Furrow1 = Ratio of average intensities within boxes of
 %%%%% Furrow1/Furrow2
 %%%%% Furrow1_S1 and Furrow1_S2 = Average intensities within the analysis
 %%%%% box
   [Furrow1,Furrow1_S1,Furrow1_S2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

    x = xFurrow2{lauf}
    y = yFurrow2{lauf}
 
 [BW_Box BWHor BWVer,imgMergeSaveVertical]=doDrawBoxAtCoordinates_Hor_1_1(x,y,imgMergeSaveVertical,voxelX_mum);
 
       
       [Furrow2,Furrow2_S1,Furrow2_S2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

    x = xPole1{lauf}
    y = yPole{lauf}
 
 [BW_Box BWHor BWVer,imgMergeSaveVertical]=doDrawBoxAtCoordinates_Ver_1_1(x,y,imgMergeSaveVertical,voxelX_mum);
 
       [Pole1, Pole1_S1, Pole1_S2] = doCalculateRatioBoundingBox(ratio,BWVer, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

 
    x = xPole2{lauf}
    y = yPole2{lauf}
 
 [BW_Box BWHor BWVer,imgMergeSaveVertical]=doDrawBoxAtCoordinates_Ver_1_1(x,y,imgMergeSaveVertical,voxelX_mum);
 
 [Pole2, Pole2_S1, Pole2_S2] = doCalculateRatioBoundingBox(ratio,BWVer, BW_Box,S1(:,:,lauf), S2(:,:,lauf))

  saveImg(:,:,:,lauf) = imgMergeSaveVertical;
    
   
     Furrow1Raw_S1_Store(lauf) = Furrow1_S1;
     Furrow1Raw_S2_Store(lauf) = Furrow1_S2;


     Furrow2Raw_S1_Store(lauf) = Furrow2_S1;
     Furrow2Raw_S2_Store(lauf) = Furrow2_S2;


     Pole1Raw_S1_Store(lauf) = Pole1_S1
     Pole1Raw_S2_Store(lauf) = Pole1_S2

     Pole2Raw_S1_Store(lauf) = Pole2_S1
     Pole2Raw_S2_Store(lauf) = Pole2_S2


     ratioFurrow = (Furrow1+Furrow2)/2;

        ratioPole = (Pole1+Pole2)/2;

        ratioFurrowStore(lauf) = ratioFurrow
        ratioPoleStore(lauf) = ratioPole
    
 
end

tiffwrite_RBG(saveImg,[tifFilename,'_ROI.tif']);

%%%%%% This section normalizes Data to pre-Anaphase
     meanFurrow1_S1_preAna = mean(Furrow1Raw_S1_Store(anaOnset-2:anaOnset));
     meanFurrow1_S2_preAna = mean(Furrow1Raw_S2_Store(anaOnset-2:anaOnset));
    
     meanFurrow2_S1_preAna = mean(Furrow2Raw_S1_Store(anaOnset-2:anaOnset));
     meanFurrow2_S2_preAna = mean(Furrow2Raw_S2_Store(anaOnset-2:anaOnset));
     
     
     meanPole1_S1_preAna = mean(Pole1Raw_S1_Store(anaOnset-2:anaOnset));
     meanPole1_S2_preAna = mean(Pole1Raw_S2_Store(anaOnset-2:anaOnset));
     
     meanPole2_S1_preAna = mean(Pole2Raw_S1_Store(anaOnset-2:anaOnset));
     meanPole2_S2_preAna = mean(Pole2Raw_S2_Store(anaOnset-2:anaOnset));
    

   
    normFurrow1_S1 = Furrow1Raw_S1_Store ./ meanFurrow1_S1_preAna;
    normFurrow1_S2 = Furrow1Raw_S2_Store ./ meanFurrow1_S2_preAna;
    
    
    normFurrow2_S1 = Furrow2Raw_S1_Store ./ meanFurrow2_S1_preAna;
    normFurrow2_S2 = Furrow2Raw_S2_Store ./ meanFurrow2_S2_preAna;
    
    normPole1_S1 = Pole1Raw_S1_Store ./ meanPole1_S1_preAna;
    normPole1_S2 = Pole1Raw_S2_Store ./ meanPole1_S2_preAna;
    
    normPole2_S1 = Pole2Raw_S1_Store ./ meanPole2_S1_preAna;
    normPole2_S2 = Pole2Raw_S2_Store ./ meanPole2_S2_preAna;
   
   
    
    Furrow_S1_Avg = (normFurrow1_S1 + normFurrow2_S1)/2;
    Furrow_S2_Avg = (normFurrow1_S2 + normFurrow2_S2)/2;
    
    Pole_S1_Avg = (normPole1_S1 + normPole2_S1)/2;
    Pole_S2_Avg = (normPole1_S2 + normPole2_S2)/2;
    
    
    %%%% Plot Data
    
    
anaTime = anaOnset*timeInterval;
timeMax = (p*timeInterval+p) - anaTime;
anaTime = 0;
 
 timeVec = 1:timeInterval:p*timeInterval;
            
timeVec = timeVec-timeVec(anaOnset);
    
ingressionTime = timeVec(ingressionFrame);
             
    
 h=figure

 plot(timeVec,  ratioFurrowStore,'-r')
 hold on
  plot(timeVec,  ratioPoleStore,'-b')

axis([timeVec(1) 500 0.7 1.6])   

yL = get(gca,'YLim');
line([ingressionTime ingressionTime],yL,'Color','r');


legend('ratio furrow', ...
  'ratio pole')


xlabel ('Time [s]','FontSize', 16);
ylabel('Ratio [A.U.]','FontSize', 16);
title(['Ratio green/red furrow and pole' tifFilename],'FontSize', 16);

print(h,'-dpdf', [tifFilename,'_Box_05_1_ratioFurrowPole.pdf']);%tifCurvetifFilename);


close all

%%%%%% Plot raw intensities


h = figure


 plot(timeVec, Furrow_S1_Avg,'-r')
 hold on
  plot(timeVec,Pole_S1_Avg,'-b')

axis([timeVec(1) 500 0 2.5])   

  
yL = get(gca,'YLim');
line([ingressionTime ingressionTime],yL,'Color','r');


legend('raw furrow', ...
  'raw pole')


xlabel ('Time [s]','FontSize', 16);
ylabel('Normalized inensities [A.U.]','FontSize', 16);
title(['Normalized to Anaphase onset furrow and pole' tifFilename],'FontSize', 16);

print(h,'-dpdf', [tifFilename,'_Box_05_NormFurrowPole.pdf']);%tifCurvetifFilename);


close all



for lauf = 1:p

    contractileRingDistanceTmp(lauf) = contractileRingDistance{lauf};
   

end




      saveVariables = {};

            saveVariables{1} = timeVec';%time{1};
            saveVariables{2} = contractileRingDistanceTmp';
            
            saveVariables{3} = ratioFurrowStore';
            saveVariables{4} = ratioPoleStore';
            saveVariables{5} = Furrow_S1_Avg';
            saveVariables{6} = Pole_S1_Avg';
           
            saveVariables{7} = Furrow1Raw_S1_Store';
            saveVariables{8} = Pole1Raw_S1_Store';
            
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, contractile_ring_distance, ratio_furrow, ratio_pole',...
                'furrow Norm, poleNorm, Furrow Raw, Pole Raw'];
            
            outid = fopen([tifFilename,'Box_RatioAnalysis_mid.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'Box_RatioAnalysis_mid.csv'],csvData,'roffset',1,'-append')

          % save([tifFilename,'.mat'])
            tifFilename
   

%clear all