%%%%%%
%%%%%% Actin measurement and normalization to anaphase onset for
%%%%%% blebbistatin treated cells. The mat file has to be loaded into the
%%%%%% workspace
%%%%%%

ratio = greenImg;

S1 = greenImg;
S2 = redImg;

for lauf = 1:AnalysisEnd
    
    S1(lauf) = S1(lauf) - greenBkg(lauf);
    S2(lauf) = S2(lauf) - redBkg(lauf);

end

S2 = S2 .*correctionFactor;


zSectionToAnalyze = 1

imgMergeCrop = imgMerge(:,:,:,cropStart:cropEnd);
clear imgMerge
imgMerge = imgMergeCrop;

for lauf = 1:AnalysisEnd
        
        analysisFrame = zSectionToAnalyze;

        [BW_Box_Out BWHor_Out BWVer_Out x_Out y_Out] = doDrawAnalysisBox12x1_ver(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        
        BW_Box = BW_Box_Out{1};
        BWHor = BWHor_Out{1};
        BWVer  =BWVer_Out{1};
        
        x = x_Out{1};
        y = y_Out{1};
        
        xFurrow1{lauf} = x;
        yFurrow1{lauf} = y;
        
        
        Furrow1_BW_BoxStore{lauf} = BW_Box;
        Furrow1_BWHor_BoxStore{lauf} = BWHor;
        Furrow1_BWVer_BoxStore{lauf} = BWVer;
        
        
        [Furrow1,Furrow1_S1,Furrow1_S2] = doCalculateRatioBoundingBox(ratio,BWVer, BW_Box, S1(:,:,lauf), S2(:,:,lauf))


       % [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        
        BW_Box = BW_Box_Out{2};
        BWHor = BWHor_Out{2};
        BWVer  =BWVer_Out{2};
        
        x = x_Out{2};
        y = y_Out{2};
       
       [Furrow2,Furrow2_S1,Furrow2_S2] = doCalculateRatioBoundingBox(ratio,BWVer, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

        xFurrow2{lauf} = x;
        yFurrow2{lauf} = y;

        Furrow2_BW_BoxStore{lauf} = BW_Box;
        Furrow2_BWHor_BoxStore{lauf} = BWHor;
        Furrow2_BWVer_BoxStore{lauf} = BWVer;
        

       % [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
       
       
        BW_Box = BW_Box_Out{3};
        BWHor = BWHor_Out{3};
        BWVer  =BWVer_Out{3};
        
        x = x_Out{3};
        y = y_Out{3};
       
       xPole1{lauf} = x;
       yPole{lauf} = y;
       
       [Pole1, Pole1_S1, Pole1_S2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

        

        Pole1_BW_BoxStore{lauf} = BW_Box;
        Pole1_BWHor_BoxStore{lauf} = BWHor;
        Pole1_BWVer_BoxStore{lauf} = BWVer;
        
        
         BW_Box = BW_Box_Out{4};
        BWHor = BWHor_Out{4};
        BWVer  =BWVer_Out{4};
        
        x = x_Out{4};
        y = y_Out{4};

    %    [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        [Pole2, Pole2_S1, Pole2_S2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box,S1(:,:,lauf), S2(:,:,lauf))

        xPole2{lauf} = x;
        yPole2{lauf} = y;

        Pole2_BW_BoxStore{lauf} = BW_Box;
        Pole2_BWHor_BoxStore{lauf} = BWHor;
        Pole2_BWVer_BoxStore{lauf} = BWVer;
    
        Furrow1Store{lauf} = Furrow1
        Furrow2Store{lauf} = Furrow2
        FurrowMeanIntensity{lauf} = (Furrow1+Furrow2)/2;

        Furrow1Raw_S1_Store(lauf) = Furrow1_S1;
        Furrow1Raw_S2_Store(lauf) = Furrow1_S2;


        Furrow2Raw_S1_Store(lauf) = Furrow2_S1;
        Furrow2Raw_S2_Store(lauf) = Furrow2_S2;

        FurrowRawMeanStore(lauf) = (Furrow2_S1 + Furrow1_S1)/2
        
         Furrow_S1_RawMeanStore(lauf) = (Furrow2_S1 + Furrow1_S1)/2
         Furrow_S2_RawMeanStore(lauf) = (Furrow2_S2 + Furrow1_S2)/2

        Pole1Raw_S1_Store(lauf) = Pole1_S1
        Pole1Raw_S2_Store(lauf) = Pole1_S2

        Pole2Raw_S1_Store(lauf) = Pole2_S1
        Pole2Raw_S2_Store(lauf) = Pole2_S2

        PoleRawMeanStore(lauf) = (Pole1_S1 + Pole1_S2)/2
         
        Pole_S1_RawMeanStore(lauf) = (Pole1_S1 + Pole2_S1)/2
        Pole_S2_RawMeanStore(lauf) = (Pole1_S2 + Pole2_S2)/2


        
        Pole1Store{lauf} = Pole1
        Pole2Store{lauf} = Pole2
        PoleMeanIntensity{lauf}= (Pole1+Pole2)/2;

        ratioFurrow = (Furrow1+Furrow2)/2;

        ratioPole = (Pole1+Pole2)/2;

        ratioFurrowStore{lauf} = ratioFurrow
        ratioPoleStore{lauf} = ratioPole
        
        
        contractileRingDistance{lauf} = (pdist2([xFurrow1{lauf},yFurrow1{lauf}],[xFurrow2{lauf},yFurrow2{lauf}])) * voxelX_mum


        
end





redMax1_New = redMax_1_new
redMax2_New = redMax_2_new

clear expV_F1_SammelTmpPlot
clear expV_F2_SammelTmpPlot




for lauf =1:p
    
       

    yGreen1BkgTmp =  yGreen1Bkg{lauf};
    yGreen2BkgTmp =  yGreen2Bkg{lauf};
    
    
    yRed1BkgTmp =  yRed1Bkg{lauf};
    yRed2BkgTmp =  yRed2Bkg{lauf};
    
    
    
        
        
         try yGreen1BkgTmpStore(lauf) = mean(yGreen1BkgTmp(redMax1_New(lauf)-10:redMax1_New(lauf)+10));
            catch ME
                warning(['Shift rmark 1.', num2str(lauf)]);
              yGreen1BkgTmpStore(lauf) = mean(yGreen2BkgTmp(redMax2_New(lauf)-10:redMax2_New(lauf)+10));
         end
        
        try yGreen2BkgTmpStore(lauf) = mean(yGreen2BkgTmp(redMax2_New(lauf)-10:redMax2_New(lauf)+10));
        catch ME
                warning(['Shift rmark 2.', num2str(lauf)]);
              yGreen2BkgTmpStore(lauf) = mean(yGreen1BkgTmp(redMax1_New(lauf)-10:redMax1_New(lauf)+10));
        end
         
        
        
        
         try yRed1BkgTmpStore(lauf) = mean(yRed1BkgTmp(redMax1_New(lauf)-10:redMax1_New(lauf)+10));
            catch ME
                warning(['Shift rmark 1.', num2str(lauf)]);
              yRed1BkgTmpStore(lauf) = mean(yRed2BkgTmp(redMax2_New(lauf)-10:redMax2_New(lauf)+10));
         end
        
        try yRed2BkgTmpStore(lauf) = mean(yRed2BkgTmp(redMax2_New(lauf)-10:redMax2_New(lauf)+10));
        catch ME
                warning(['Shift rmark 2.', num2str(lauf)]);
              yRed2BkgTmpStore(lauf) = mean(yRed1BkgTmp(redMax1_New(lauf)-10:redMax1_New(lauf)+10));
         end
        
        
        
end




  
furrowRatioSpline = ((yGreen1BkgTmpStore ./ (yRed1BkgTmpStore .* correctionFactor)) +  (yGreen2BkgTmpStore ./ (yRed2BkgTmpStore .* correctionFactor)))/2;

intenistyFurrowSpline_A = (yGreen1BkgTmpStore + (yRed1BkgTmpStore .* correctionFactor)) ./ 2
intenistyFurrowSpline_B = (yGreen2BkgTmpStore + (yRed2BkgTmpStore .* correctionFactor)) ./ 2

intensityFurrowSpline = (intenistyFurrowSpline_A + intenistyFurrowSpline_B) ./ 2


furrowRatio = (Furrow_S1_RawMeanStore)./(Furrow_S2_RawMeanStore);
poleRatio = Pole1Raw_S1_Store ./ Pole1Raw_S2_Store;


intensityFurrow = (Furrow_S1_RawMeanStore + Furrow_S2_RawMeanStore) ./ 2;
intensityPole = (Pole_S1_RawMeanStore + Pole_S2_RawMeanStore) ./ 2

%%%%%
%%%%% normalize intensities to anaphase onset
%%%%%
anaOnset = find(0 == timeVec);

intensityFurrowSplineCrop = intensityFurrowSpline(2:end);



if anaOnset == 1
    
    
    intensityFurrowNorm = intensityFurrow ./ mean(intensityFurrow(anaOnset))
    intensityPoleNorm = intensityPole ./  mean(intensityPole(anaOnset))
    intensityFurrowSplineNorm  = intensityFurrowSpline ./ mean(intensityFurrowSpline(anaOnset))
    intensityFurrowSplineNormCrop  = intensityFurrowSplineCrop ./ mean(intensityFurrowSplineCrop(anaOnset))

elseif anaOnset == 2
    
    intensityFurrowNorm = intensityFurrow ./ mean(intensityFurrow(anaOnset-1:anaOnset))
    intensityPoleNorm = intensityPole ./  mean(intensityPole(anaOnset-1:anaOnset))
    intensityFurrowSplineNorm  = intensityFurrowSpline ./ mean(intensityFurrowSpline(anaOnset-1:anaOnset))
    intensityFurrowSplineNormCrop  = intensityFurrowSplineCrop ./ mean(intensityFurrowSplineCrop(anaOnset-1:anaOnset))

elseif anaOnset > 2
    
    intensityFurrowNorm = intensityFurrow ./ mean(intensityFurrow(anaOnset-2:anaOnset))
    intensityPoleNorm = intensityPole ./  mean(intensityPole(anaOnset-2:anaOnset))
    intensityFurrowSplineNorm  = intensityFurrowSpline ./ mean(intensityFurrowSpline(anaOnset-2:anaOnset))
    intensityFurrowSplineNormCrop  = intensityFurrowSplineCrop ./ mean(intensityFurrowSplineCrop(anaOnset-1:anaOnset))

end
close all

plot(timeVec, intensityFurrowNorm(1:p),'r')
hold on
plot(timeVec, intensityPoleNorm(1:p),'b')
plot(timeVec, intensityFurrowSplineNorm(1:p),'k')
%plot(timeVec(2:end), intensityFurrowSplineNormCrop(1:p),'k')




    saveVariables = {};

            saveVariables{1} = timeVec(1:end)';
            saveVariables{2} =  intensityFurrowSplineNorm(1:p)';
            
            saveVariables{3} =   intensityPoleNorm(1:p)';
              
           
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, intensityFurrow, intensityPole ']
            
            outid = fopen([folderName,'160802_actinIntensity.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([folderName,'160802_actinIntensity.csv'],csvData,'roffset',1,'-append')
           
    varOut = [timeVec' intensityFurrowSplineNorm(1:p)' intensityPoleNorm(1:p)']
    
    open varOut
 %%               
    save([folderName,'intensities.mat'])
    
      clear all
            close all


