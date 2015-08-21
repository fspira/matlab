%%%%%% Generate nice anisotropy colormap for 5 frame poslcope data - load
%%%%%% workspace of the respective cell



segment = doSegmentImage( greenS);


for lauf = 1:p
    L2 = segment(:,:,lauf);
    cc = regionprops(L2,'all');
    B = bwboundaries(L2,8);
    boundary = B{1};
    
    bStore{lauf} = B{1};
    
    MIJ.createImage(I0File(:,:,lauf));
    MIJ.run('setLine8');
    MIJ.setRoi( [ boundary(:,2)'-112; boundary(:,1)'-112], ij.gui.Roi.POLYLINE);
  
    k = waitforbuttonpress
    clear logData
    
    MIJ.run('getSplineCoordinates')
    logData = MIJ.getLog;
    MIJ.run('closeLogWindow')
    logDataTmp = char(logData);
    logDataTmp1 = strread(logDataTmp, '%s');
    
    xLogTmp = logDataTmp1(2:4:end);
    yLogTmp = logDataTmp1(3:4:end);
    
    
    xLog = zeros(length(xLogTmp),1,'double');
    yLog = zeros(length(yLogTmp),1,'double');
    
    for subrun = 1 : length(xLogTmp)
        
        xLog(subrun) = round(str2num(xLogTmp{subrun}));
        yLog(subrun) = round(str2num(yLogTmp{subrun}));
        
    end
    
    newLinescan{lauf} = [xLog yLog]
    
    
    
    
    MIJ.run('closeAllWindows');
end


            
                
for lauf = 1:AnalysisEnd %length(redStack)
        
    
    
    
     
    %%%% Load XY coordinates of the two flanking regions into the tmp
    %%%% variable which is later used for analysis. This variable will be
    %%%% updated every loop
         
         flank1Tmp = newLinescan{lauf};
        
         
         %%%% Load the red image
         MIJ.createImage(I0File(:,:,lauf));
         MIJ.run('setLine12');%%%% increase linescan width
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);%%%% Draw line into the image
         
         
         %%%%% Update InterpolatedLine
          clear logData
          MIJ.run('getSplineCoordinates')%%%% Get interpolated coordinates of the line
          logData = MIJ.getLog;%%%%% Reading MIJI log file to get coordinates
          
          MIJ.run('closeLogWindow')%%%% Close log file
         
          MIJ.run('getLinescanRed');%%%% Get intensities under the line
         yRed1{lauf} = MIJ.getColumn('y');%%%% Intensity values
         xRed1{lauf} = MIJ.getColumn('x');%%%% Running variable for each pixel can be discarded
        logDataTmp = char(logData);%%%%% convert values to characters - important to search for a string in the next step
        logDataTmp1 = strread(logDataTmp, '%s'); %%%% Sort log data using regular expression
        
         xLogTmp = logDataTmp1(2:4:end); %%%% Save X coordinates of the interpolated linescan
         yLogTmp = logDataTmp1(3:4:end);%%% Save Y coordinates of the interpolatd linescan
         
             
        xLog = zeros(length(xLogTmp),1,'double'); %%%% Pre allocate X and Y variables
        yLog = zeros(length(yLogTmp),1,'double');
         
         for subrun = 1 : length(xLogTmp)%%%% Save XY coordinates of the interpolated linescan
             
             xLog(subrun) = round(str2num(xLogTmp{subrun}));
             yLog(subrun) = round(str2num(yLogTmp{subrun}));
             
         end
         
         flank1StoreNew{lauf} = [xLog yLog]%%% Save coordinates into store array
       
        
             MIJ.run('closeAllWindows');
           
            % flank1Tmp = flank1StoreNew{lauf};
         
             %%%%%% Get lineprofiles from the green channel
         MIJ.createImage(I0File(:,:,lauf));
         MIJ.run('setLine12');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI0_1{lauf} = MIJ.getColumn('y');
         xI0_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
              MIJ.createImage(I45File(:,:,lauf));
         MIJ.run('setLine12');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI45_1{lauf} = MIJ.getColumn('y');
         xI45_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
             MIJ.createImage(I90File(:,:,lauf));
         MIJ.run('setLine12');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI90_1{lauf} = MIJ.getColumn('y');
         xI90_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
              MIJ.createImage(I135File(:,:,lauf));
         MIJ.run('setLine12');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI135_1{lauf} = MIJ.getColumn('y');
         xI135_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
          
              MIJ.createImage(IBFile(:,:,lauf));
         MIJ.run('setLine12');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yIB_1{lauf} = MIJ.getColumn('y');
         xIB_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
        
         
        
         
         
end



slidingWindow =7




for lauf = 1:p
   
    [anisoOut] = doGenerateSlidingAnsiotropy(yI0_1{lauf},yI45_1{lauf},yI90_1{lauf},yI135_1{lauf},slidingWindow,ItoSMatrix)
    
    anisoOut_Store1{lauf} = anisoOut
    
    
end


cmap = jet;
 cmap(1,3) = 0;

for lauf = 1:p




file1Tmp = I0File(:,:,lauf);
file1Tmp(:,:) = 0;
file1Tmp = double(file1Tmp);
file2Tmp = file1Tmp;

flank1Tmp = newLinescan{lauf};

cFurrowContour1 =  flank1Tmp(:,1);
rFurrowContour1 =  flank1Tmp(:,2);


cFurrowContour1 = doEliminateNegativeValues(cFurrowContour1)
rFurrowContour1 = doEliminateNegativeValues(rFurrowContour1)
cFurrowContour1 = minDistBoundary(cFurrowContour1)
rFurrowContour1 = minDistBoundary(rFurrowContour1)


anisoOutPlot1 = anisoOut_Store1{lauf};

for i = 1:length(rFurrowContour1)
   
    file1Tmp(rFurrowContour1(i),cFurrowContour1(i)) = anisoOutPlot1(i);
    
end


file1Norm = normalizedImage(file1Tmp);

imshow(file1Norm,[])


 
 
 [res1,res_tmp1]=splineDilate(file1Tmp,file1Norm,[cFurrowContour1 rFurrowContour1], [1,3], 3);
 

 
 
resNorm = normalizedImage(res1);
resStoreNorm(:,:,lauf) = resNorm;
resStore(:,:,lauf) = res1;
 
imshow(resNorm,[])
colormap(cmap)



pause(0.5)

end

curdir = pwd;

tiffwrite_mat(resStoreNorm,[curdir,'/analysis/anisotropyContourNorm_SlidingWindow7_lineWidth12'])


tiffwrite_mat(resStore,[curdir,'/analysis/anisotropyContourOrig_SlidingWindow7_lineWidth12'])

 save(['workspace','_Linescan.mat'])
