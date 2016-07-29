

%%%%%% Recompute bkg
clear img


%javaaddpath  'C:\Programme\fiji-win64-20110307\Fiji.app\scripts\Miji.m'
%javaaddpath 'C:\Programme\fiji-win64-20110307\Fiji.app\jars\ij-1.47d.jar'
%addpath  'C:\Programme\fiji-win64-20110307\Fiji.app\scripts\'
%addpath 'C:\Programme\fiji-win64-20110307\Fiji.app\macros'
%addpath('C:\Users\spira\Dropbox\MatlabSkripte\Image_Processing_utils')
%addpath('C:\Users\spira\Dropbox\MatlabSkripte\SubPixel')




javaaddpath '/Applications/Fiji.app/scripts/Miji.m'
javaaddpath '/Applications/Fiji.app/jars/ij-1.48s.jar'

javaaddpath '/Applications/MATLAB_R2012a_Student.app/java/mij.jar'


addpath('/Applications/Fiji.app/scripts')

addpath('/Applications/Fiji.app/plugins')
addpath('/Users/spira/Documents/Matlab_scripte/Image_Processing_utils')

addpath('/Users/spira/Documents/Matlab_scripte/tiffIO')
addpath('/Users/spira/Documents/Matlab_scripte/')
addpath('/Users/spira/Desktop/Desktop/LifeactCherry_GlGPIEgfp/131204')
addpath('/Users/spira/Documents/MATLAB_scripte/ImageProcessing/Utilities')

addpath('/Users/spira/Desktop/programme/calculateAnisotropieRatio')
addpath('/Users/spira/Desktop/programme/centerOfMass')
addpath('/Users/spira/Desktop/programme/curvature')
addpath('/Users/spira/Desktop/programme/determineAngle')
addpath('/Users/spira/Desktop/programme/staging')
addpath('/Users/spira/Desktop/programme/tools')

addpath('/Users/spira/Desktop/programme/mijiTools')

curdir = pwd;

%MIJ.setupExt ('/Users/spira/Desktop/Desktop/LifeactCherry_GlGPIEgfp/131204');

   
     
  

%bleachCorrGreen = load('/Users/spira/Desktop/02330/131126/ActinGFP_MyrPalmCherry/bleachCorr/ActinGFP_MyrPalmCherry_stack1-1.tif_BleachCorrectionActinGFP.mat');

%bleachCorrGreen = bleachCorrGreen.yGreen;



%bleachCorrRed =  load('/Users/spira/Desktop/02330/131126/ActinGFP_MyrPalmCherry/bleachCorr/ActinGFP_MyrPalmCherry_stack1-2.tif_BleachCorrectionMyrPalmCherry.mat');


%bleachCorrRed = bleachCorrRed.yRed;

clear redStack
clear greenStack

imgOrigGreen = tiffread30([char(tifFilename(1:truncName-1)), '_greenCrop.tif']);
imgOrigRed = tiffread30([char(tifFilename(1:truncName-1)), '_redCrop.tif']);




greenStack = cat(3,imgOrigGreen.data);
redStack = cat(3,imgOrigRed.data);


%tiffwrite_mat(greenStack, [folderName,'_greenCrop.tif']);
%tiffwrite_mat(redStack, [folderName,'_redCrop.tif']);

cd([curdir,'/','averagedBackground'])

anaOnset =11;


clear imgOrigGreen;
clear imgOrigRed;
 
[m n p] = size(greenStack);

%%%%%%%%%%%%%%% background correction and smoothing

[bkga bkgb bkg1 bkg2 bkg3 bkg4,imgOut] = bkgCorrectionRedGreenOutsideInside_manual(greenStack, redStack,voxelX_mum);

%%%%%%%%%%%%%% bleach correction

%bleachCorrGreen = bleachCorrGreen(1:p);
%bleachCorrRed = bleachCorrRed(1:p);

greenStack = double(greenStack);
redStack = double(redStack);
curdir = pwd;
Miji;
cd(curdir)
%AnalysisEnd = p;

for lauf = 1:AnalysisEnd %length(redStack)
        
       
         
         
         MIJ.createImage(redStack(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ cNew1Ref{lauf}'; rNew1Ref{lauf}'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRed1{lauf} = MIJ.getColumn('y');
         xRed1{lauf} = MIJ.getColumn('x');
        MIJ.run('closeAllWindows');
         
         MIJ.createImage(redStack(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ cNew2Ref{lauf}'; rNew2Ref{lauf}'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRed2{lauf} = MIJ.getColumn('y');
         xRed2{lauf} = MIJ.getColumn('x');
        MIJ.run('closeAllWindows');
         
         
         
         
         MIJ.createImage(redStack(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ cNonPol1{lauf}'; rNonPol1{lauf}'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedNonPol1{lauf} = MIJ.getColumn('y');
         xRedNonPol1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         MIJ.createImage(redStack(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ cNonPol2{lauf}'; rNonPol2{lauf}'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedNonPol2{lauf} = MIJ.getColumn('y');
         xRedNonPol2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
          
         
         
         MIJ.createImage(greenStack(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ cNew1Ref{lauf}'; rNew1Ref{lauf}'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yGreen1{lauf} = MIJ.getColumn('y');
         xGreen1{lauf} = MIJ.getColumn('x');
        MIJ.run('closeAllWindows');
         
         MIJ.createImage(greenStack(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ cNew2Ref{lauf}'; rNew2Ref{lauf}'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yGreen2{lauf} = MIJ.getColumn('y');
         xGreen2{lauf} = MIJ.getColumn('x');
        MIJ.run('closeAllWindows');
         
         
         
         
         MIJ.createImage(greenStack(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ cNonPol1{lauf}'; rNonPol1{lauf}'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yGreenNonPol1{lauf} = MIJ.getColumn('y');
         xGreenNonPol1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         MIJ.createImage(greenStack(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ cNonPol2{lauf}'; rNonPol2{lauf}'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yGreenNonPol2{lauf} = MIJ.getColumn('y');
         xGreenNonPol2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
         
end
%p=42;
%%%%%%% Background correct Values
yRed1_orig =yRed1;
yRed2_orig = yRed2;

yGreen1_orig =yGreen1;
yGreen2_orig = yGreen2;

yGreenNonPol1_orig = yGreenNonPol1;
yGreenNonPol_orig = yGreenNonPol2;


yRedNonPol1_orig = yRedNonPol1;
yRedNonPol2_orig = yRedNonPol1;


for lauf =1:p
    yRed1{lauf} = yRed1{lauf}-double(bkgb(lauf));
    yRed2{lauf} = yRed2{lauf}-bkgb(lauf);

    yGreen1{lauf}  =yGreen1{lauf} -bkga(lauf);
    yGreen2{lauf} = yGreen2{lauf} - bkga(lauf);

    yRedNonPol1{lauf} = yRedNonPol1{lauf} - bkgb(lauf);
    yRedNonPol12{lauf}= yRedNonPol2{lauf} -bkgb(lauf)


    yGreenNonPol1{lauf} = yGreenNonPol1{lauf} -bkga(lauf);
    yGreenNonPol2{lauf} = yGreenNonPol2{lauf} -bkga(lauf);
end

%test = yRedMid1

%%%%%%% Identification of midsections of the ROI's. The short script plots
%%%%%%% a high pixel value into the linescan, which can be easily detected
%%%%%%% as a peak and used to find the mid regin.

%Miji;


%%%%%% All intensities calculated for left and right flank


%%%% set sliding window size in micron
windowSize = 4; %%% in microns

sWSet = round(windowSize / voxelX_mum);


%sWSet = 40;
sW = sWSet;
sW1 = sW; %%% in pixels
sW2 = sW;
clear  greenValuesSliding_1
clear  greenValuesSliding_2
clear  redValuesSliding_1
clear  redValuesSliding_2

for lauf = 1:AnalysisEnd    

    redValuesSliding_1tmp = yRed1{lauf};
    
        %%%%% This section ensures that the sliding windows fits the
        %%%%% selection. If not the selection is cropped to fit the window
    
    %    redMax1(dist1 == 1) = 2;
        
        if redMax1(lauf)+ (sW1/2) > length(redValuesSliding_1tmp)
            
            sW1 = round(length(redValuesSliding_1tmp)-(redMax1(lauf)))-1
            if  mod(sW1,2) ==1
                sW1 = sW1 -1;
            end
            
        elseif redMax1(lauf)-(sW1/2) < 1
            sW1 = round(redMax1(lauf)/2)-1
            if  mod(sW1,2) ==1
                sW1 = sW1 -1;
            end
        end
        
    redValuesSliding_1{lauf} = redValuesSliding_1tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
    

    redValuesSliding_2tmp = yRed2{lauf};
    
     if redMax2(lauf)+ (sW2/2) > length(redValuesSliding_2tmp)
            
            sW2 = round(length(redValuesSliding_2tmp)-(redMax2(lauf)))-1
            if  mod(sW2,2) ==1
                sW2 = sW2 -1;
            end
        elseif redMax2(lauf)-(sW2/2) < 1
             sW2 = round(redMax2(lauf)/2)-1
             if  mod(sW2,2) ==1
                sW2 = sW2 -1;
            end
        end
    sWStore1(lauf) = sW1
      sWStore2(lauf) = sW2
    

    redValuesSliding_1tmp = yRed1{lauf};
    redValuesSliding_1{lauf} = redValuesSliding_1tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   
    
    redValuesSliding_2tmp = yRed2{lauf};
    redValuesSliding_2{lauf} =   redValuesSliding_2tmp([redMax2(lauf)-(sW2/2):redMax2(lauf)+(sW2/2)]);
   
    
    
    
    greenValuesSliding_1tmp = yGreen1{lauf};
    greenValuesSliding_1{lauf} = greenValuesSliding_1tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   
    
    greenValuesSliding_2tmp = yGreen2{lauf};
    greenValuesSliding_2{lauf} =   greenValuesSliding_2tmp([redMax2(lauf)-(sW2/2):redMax2(lauf)+(sW2/2)]);
   
    
    
    
    %%%%% Window which selects 80 pixel of the pole Roi.
    %%%% Select Window Size
    
    SW = sW;
    SW = SW/2; 
    
    
     
    %%%%%%%%%% Merge pole Regions red
    
     tmpyRedNonPol1 =    yRedNonPol1{lauf};
  
     
    yRedNonPol1Cut{lauf} =    tmpyRedNonPol1([redPoleMax1(lauf)-(SW/2):redPoleMax1(lauf)+(SW/2)]);
  
     
    tmpyRedNonPol2 =    yRedNonPol2{lauf};
    
      yRedNonPol2Cut{lauf} =    tmpyRedNonPol2([redPoleMax2(lauf)-(SW/2):redPoleMax2(lauf)+(SW/2)]);
    

    redValuesNonIngMean(lauf) = (mean(yRedNonPol1Cut{lauf}) + mean(yRedNonPol2Cut{lauf}))/2;
    
    %%%%%%%%%% Merge pole regions green
      
    tmpyGreenNonPol1 =   yGreenNonPol1{lauf};
  
     
    yGreenNonPol1Cut{lauf} =   tmpyGreenNonPol1([redPoleMax1(lauf)-(SW/2):redPoleMax1(lauf)+(SW/2)]);
  
     
    tmpyGreenNonPol2 =    yGreenNonPol2{lauf};
    
      yGreenNonPol2Cut{lauf} =   tmpyGreenNonPol2([redPoleMax2(lauf)-(SW/2):redPoleMax2(lauf)+(SW/2)]);
    

   greenValuesNonIngMean(lauf) = (mean(yGreenNonPol1Cut{lauf}) + mean(yGreenNonPol2Cut{lauf}))/2;
    
    
    
    %%%%%%%%%%%%%%%%%%
    
    
   
   
    %%%%%% Mean intensity per sliding windows
  
    redValuesSlidingMean_1(lauf) = mean(redValuesSliding_1{lauf});
    redValuesSlidingMean_2(lauf) = mean(redValuesSliding_2{lauf});
    
    
    
    greenValuesSlidingMean_1(lauf) = mean(greenValuesSliding_1{lauf});
    greenValuesSlidingMean_2(lauf) = mean(greenValuesSliding_2{lauf});
    
    
    %%%%%%% Average both flanks
  
     
     redValuesSlidingFlankMean(lauf) = (mean(redValuesSliding_1{lauf}) + mean(redValuesSliding_2{lauf}))/2;
    
     greenValuesSlidingFlankMean(lauf) = (mean(greenValuesSliding_1{lauf}) + mean(greenValuesSliding_2{lauf}))/2;
    
     
    sW = sWSet;
sW1 = sW; %%% in pixels
sW2 = sW;

end

%%%%%% Normalize linescans to anaphase onset


    if anaOnset  >= 3

          redFurrowNorm = redValuesSlidingFlankMean ./mean(redValuesSlidingFlankMean(anaOnset-2:anaOnset));
          greenFurrowNorm = greenValuesSlidingFlankMean ./mean(greenValuesSlidingFlankMean(anaOnset-2:anaOnset));

         greenValuesPoleNorm = greenValuesNonIngMean ./mean(greenValuesNonIngMean(anaOnset-2:anaOnset));

         redValuesPoleNorm = redValuesNonIngMean ./mean(redValuesNonIngMean(anaOnset-2:anaOnset));

    else


       redFurrowNorm = redValuesSlidingFlankMean ./mean(redValuesSlidingFlankMean(1:anaOnset));
       greenFurrowNorm = greenValuesSlidingFlankMean ./mean(greenValuesSlidingFlankMean(1:anaOnset));

       
         greenValuesPoleNorm = greenValuesNonIngMean ./mean(greenValuesNonIngMean(1:anaOnset));
       
       
          redValuesPoleNorm = redValuesNonIngMean ./mean(redValuesNonIngMean(1:anaOnset));
       
    end



%%%%%%%% Calculate ratios



FurrowVsPoleRed = redFurrowNorm ./ redValuesPoleNorm;
FurrowVsPoleGreen = greenFurrowNorm ./ greenValuesPoleNorm;





   % cd([curdir '/' ]);

   savefile = 'ROIs.mat';
    
     

            saveVariables = {};
 
            saveVariables{1} = timeVec';
           
            saveVariables{2} = redFurrowNorm';
              saveVariables{3} = greenFurrowNorm';
              
            
           
            
            clear csvData;
            csvData=saveVariables;
   % header= ['1,2,3,4,5,6,7,8,9,10,11,12,13,14'];
    header= ['Time', 'Ratio Furrow vs. Pole Red','Ratio Furrow vs. Pole Green'];         
            
            outid = fopen('Ratios.csv', 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite('Ratios.csv',csvData,'roffset',1,'-append')
            
            
            %%%%% Save varibles from the summed flanks
  %    clear all
      

        
% load('Workspace1.mat')
 
 saveVariables = {};
          
            saveVariables{1} = timeVec';
            
            saveVariables{2} = redValuesSlidingFlankMean';
            saveVariables{3} = redValuesNonIngMean';
             saveVariables{4} = greenValuesSlidingFlankMean';
              saveVariables{5} = greenValuesNonIngMean';
          
           
          
           
            
            
            clear csvData;
            csvData=saveVariables;
   % header= ['1,2,3,4,5,6,7,8,9,10,11,12,13,14'];
            header= ['Time , redChannelFurrow, redChannelPole,greenChannelFurrow,greenChannelPole'];
            
            outid = fopen('FurrowIngressionCorticalValuesSUMFlanks.csv', 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ('FurrowIngressionCorticalValuesSUMFlanks.csv',csvData,'roffset',1,'-append')
           
            
            %%%%%% Save background
            
 saveVariables = {};
          
            saveVariables{1} = timeVec';
            
            saveVariables{2} = bkg3';
            saveVariables{3} = bkg4';
             
           
          
           
            
            
            clear csvData;
            csvData=saveVariables;
   % header= ['1,2,3,4,5,6,7,8,9,10,11,12,13,14'];
            header= ['Time , greenChannelBkg, redChannelBkg'];
            
            outid = fopen('Bkg.csv', 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ('Bkg.csv',csvData,'roffset',1,'-append')
           
            
          %%%%%% Normalize ingressing Furrow Plot
          
            save('Workspace1DoubleBackground')
            



