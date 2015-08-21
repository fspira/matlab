%%
clear all

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;

blue = zeros(256,3);
blue(:,3) = distvec;

javaaddpath '/Applications/Fiji.app/scripts/Miji.m'
javaaddpath '/Applications/Fiji.app/jars/ij-1.47b.jar'
addpath('/Applications/Fiji.app/scripts')

addpath('/Applications/Fiji.app/scripts')
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


curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell1_ana_ablation.lsm';
tifFilenameMid = 'cell1_ana_equator.lsm';



saveFileName = 'cell1_ana_save';

truncName = findstr(tifFilename,'.lsm');
emptyFlag = isempty(truncName);

if emptyFlag ==1

 truncName = findstr(tifFilename,'.tif');

end


folderName = [tifFilename(1:truncName-1),'_analysis'];
mkdir([curdir '/' folderName]);


%%%% Load image and extract metadata
imgMid  = tiffread30(tifFilenameMid);
midVoxelSize = imgMid.lsm.VoxelSizeX * 1000000;

imgMid = cat(3,imgMid.data);
imgEquator = imgMid{2};
imgMid = imgMid{1};


image = tiffread30(tifFilename);
timeInterval = getfield(image,'lsm','TimeOffset');
ablationTime = timeInterval(2);
timeInterval = timeInterval(4) - timeInterval(3);

voxelX = getfield(image,'lsm','VoxelSizeX')* 1000000;





origImage = cat(3,image.data);
[m n p] = size(origImage)

index=1
for lauf = 1:p

    origImageTmp(:,:,lauf) = origImage{index};

    index = index+2;
    
end

clear origImage;
origImage = origImageTmp;
clear origImageTmp;

%imshow(recoilImage(:,:,2),[])


cd([curdir '/' folderName]);



%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% End of header, start analysis
%%%%%%%%%%%%%%%%%

[m n p] = size(origImage)

recoilImage = origImage(:,:,2:11);

[m n p] = size(recoilImage)


%%%%%%% detect bleached area and obatin paraemters for the ellipse
close all



[recoilImageBKG bkgRoi]  = bkgCorrection(recoilImage);






%imageFilter = gradmag;
[I2 rect] = imcrop(recoilImageBKG(:,:,1),[]);
close all
%imshow(I)
%imshow(recoilImageBKG(:,:,lauf),[])
%  imshow(imageFilter(:,:,lauf))
%hold on
%    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)


%%% Set marker for seeded Watershed

bkgRoiNorm = bkgRoi .*255;

markerImage =  uint8(bkgRoiNorm)  + uint8(recoilImageBKG(:,:,1));

 imshow(markerImage(:,:,1))
    
     [xxx,yyy] = ginput(1)
     
     close all
    %% 
%     figure(1)
for lauf =1:p
    
    
       % I_crop =    imcrop(imageFilter(:,:,lauf),rect);
        
%[a]= doSegmentImage(I)
%imshow(a,[])
 I_eq(:,:,lauf) = adapthisteq(recoilImageBKG(:,:,lauf));
%   [I_eq, T] = histeq(recoilImageBKG(:,:,lauf),64) ;
   
%imshow(I_eq)
I_eq(:,:,lauf) = double(I_eq(:,:,lauf));


            h = fspecial('gaussian',[15 15],4) ;
           % imageFilter(:,:,lauf) = imfilter(recoilImageBKG(:,:,lauf),h);
             imageFilter(:,:,lauf) = imfilter(I_eq(:,:,lauf),h);

            
            mask = imageFilter(:,:,lauf);

        marker = bkgRoi;
       % marker(129:135,139:145) = true;

        J = mask;
%J(marker) = 255;
       % figure, imshow(J); title('Marker Image Superimposed on Mask');

    K = imimposemin(mask,marker);
    
    imshow(K)
    
    %  [xx,yy] = ginput(1)
    
    propsMarker = regionprops(marker,'Centroid');
    
  % xx = propsMarker.Centroid(1) - 3;
  %  yy = propsMarker.Centroid(1) - 6;
    
    
    pixelMarker = K(round(yyy),round(xxx));
    
    I2 = imhmin(K,pixelMarker);
     imshow(I2)
      % [xx,yy] = ginput(1)
    % pixelMarker = 40;
     pixelMarker = (I2(round(yyy),round(xxx)));
     
     I3= imextendedmin(I2,pixelMarker,8);
    
     regionArea = regionprops(I3,'Area')
     flag = 0;
     
         incrementValue = 1;
     
     if  regionArea.Area > 6000
     
     figure(4)
         imshow(I3)

             while regionArea.Area >1200
                 
                 flag = 0;

                 pixelMarker = (I2(round(yyy),round(xxx))) - incrementValue;


                 I3= imextendedmin(I2,pixelMarker,8);

                 regionArea = regionprops(I3,'Area');

                 incrementValue = incrementValue + 1


             end

     end
         
           incrementValue = 1;
             while regionArea.Area < 1300
                 
                 flag = 1;

                 pixelMarkerNew = pixelMarker + incrementValue;


                 I3= imextendedmin(I2,pixelMarkerNew,8);

                 regionArea = regionprops(I3,'Area');

                 incrementValue = incrementValue + 1


             end
         
             
     
            
     
    if flag  == 1
            
        pixelMarker = (pixelMarkerNew-0);


                I3= imextendedmin(I2,pixelMarker,8);

                regionArea = regionprops(I3,'Area');
                 
       
     end
%              
             
   
     
    
    
    imshow(I3)
    
    segmentedImage(:,:,lauf) = I3;
 
    
  
    
   % segmentedImage1 = I2 ==  pixelMarker;
   % segmentedImage2 = I2 ==  pixelMarker+1;
   % segmentedImage3= I2 ==  pixelMarker+1;
   %  segmentedImage4 = I2 >  pixelMarker & I2 < pixelMarker+15 ;
    
   % segmentedImage(:,:,lauf) = segmentedImage1 + segmentedImage2 + segmentedImage3 + segmentedImage4;
    
    figure(3)
    imshow(segmentedImage(:,:,lauf))
    
    %imshow(segmentedImage1)
    
    I_crop = segmentedImage(:,:,lauf);
        
   % [output BW1] = CORF(I_crop,5, 0.5);
       se = strel('disk',5);
       
       
         [BW,thresh] = edge(  I_crop,'canny');
          
       
         BW1 = edge(I_crop,'canny', thresh,6);
   % BW2 = BW1;
    BW1 = bwareaopen(BW1,15);
    
   
            BW2 = imclose(BW1,se);
           % Io = imclose(BW1, se);
       % imshow(BW2)
              
         imshow(BW1)
         
         
            

            close all

%imshow(BW1)

            %BW2 = edge(imageFilter,'canny');




            [B L]= bwboundaries(BW2,8);
            clear horizontal
            clear  boundaryLength
          imshow(BW1)
            hold on
            for k = 1:length(B)
                boundary = B{k};
            
             plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2);
                boundaryLength(k) = length(B{k});
            end


            maxBoundary = max(boundaryLength);


            maxBoundaryPosition = find(boundaryLength == maxBoundary);

            boundary = B{maxBoundaryPosition};
            %%%%% Test
            
           % se = strel('disk', 10);
          
         
            
         
           cropOriginal = I_crop;%imcrop(recoilImage(:,:,lauf),rect);
         
           
            imshow(cropOriginal)

            cropOriginal_save = cropOriginal;
            
            hold on
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2);


            blackImage = I_crop;
            blackImage(:,:) = 0;
            
            for subrun = 1:length(boundary)
                blackImage(boundary(subrun,1),boundary(subrun,2)) = 255;
                cropOriginal_save(boundary(subrun,1),boundary(subrun,2)) = 255;
            end
            
            imwrite(cropOriginal_save,['CountourOriginal_',num2str(lauf),'.tif'],'tiff');
            
            imageFilterSave(:,:,lauf) = recoilImage(:,:,lauf);
            
            
            for subrun = 1:length(boundary)
                 
             imageFilterSave(boundary(subrun,1),boundary(subrun,2),lauf) = 255;
                %cropOriginal_save(boundary(subrun,1),boundary(subrun,2)) = 255;
            end
            
          
            
            %imshow(blackImage);
            hold on
            
            ablationContour = im2bw(blackImage);
            
           % pause(0.2)
           
            
            STATS = regionprops(ablationContour, 'all');
          
            
            %  doDrawEllipse(STATS, 40,1,I_crop) 
            
            close all
            
           
            
            
            %%%%% Find points on the contour which are at 0 or 90 degree relative to ellipse orientation 
            
             close all 
            flag = 0
            
             clear angleStore
             
               orientationEllipse = STATS.Orientation ;
               
           
               
                
               contourRegularity(lauf) =   (STATS.FilledArea/STATS.Area^2);
                contourCircularity(lauf) =  STATS.EquivDiameter;
             
              
               
               if orientationEllipse > 0
               
                        flag = 1
                        
                     for subrun = 1:length(boundary)

                            dy = boundary(subrun,2) - round(STATS.Centroid(1));
                            dx = boundary(subrun,1) - round(STATS.Centroid(2));

                            plot(boundary(subrun,2),boundary(subrun,1),'xr')
                             plot(STATS.Centroid(1),STATS.Centroid(2),'xr')

                        %    dx = cChromatin1(1) - cChromatin2(1);
                         %   dy = rChromatin1(1) - rChromatin2(1);

                         
                         
                             if dy < 0

                                   angle = (atan2(dx,dy))*180/pi;

                             else

                                  angle = (atan2(dx,dy))*180/pi;

                             end


                             angleStore(subrun) = round(angle);

                     end
               
               
                else
                   
                   
                   flag = 0
                   
                    for subrun = 1:length(boundary)

                        dy = boundary(subrun,2) - round(STATS.Centroid(1));
                        dx = boundary(subrun,1) - round(STATS.Centroid(2));

                        plot(boundary(subrun,2),boundary(subrun,1),'xr')
                         plot(STATS.Centroid(1),STATS.Centroid(2),'xr')

                    %    dx = cChromatin1(1) - cChromatin2(1);
                     %   dy = rChromatin1(1) - rChromatin2(1);

                         if dy < 0

                               angle = (atan2(dy,dx))*180/pi;

                         else

                              angle = (atan2(dy,dx))*180/pi;

                         end


                         angleStore(subrun) = round(angle);


                    end
               end

            
            %%%% Map coordinates onto the contour
            
            
            
              
            
                 horizontal_A = orientationEllipse;
                
                %%%%% Find degrees which are +-90 and -180 degrees relative
                %%%%% to ellipse orientation. Ellipse orientation as an
                %%%%% output from regionprops gives -180+180 degrees. This
                %%%%% means that the desired orientations have to be mapped
                %%%%% into -180+180 system. (Avoid degrees greater -180+180
                
               % STATS.Orientation > 0
                     
                    
                    

                     if horizontal_A - 180 < - 180


                         horizontal_B = 180 + ((horizontal_A - 180) + 180);



                     else

                         horizontal_B = horizontal_A - 180;


                     end


                     if horizontal_A - 90 < - 180;


                         horizontal_C = 90 + ((horizontal_A - 180) + 180);



                     else

                         horizontal_C = horizontal_A - 90;


                     end


                     if horizontal_A + 90 >  180


                         horizontal_D = ((horizontal_A - 180) - 90);



                     else

                         horizontal_D = horizontal_A + 90;

                     
                     end
             
                 
              
           
             
                 %%%%%% Find points on the contour which are +-90 and -180
                 %%%%%% relative to the ellipse orientation. 
                 
                
            h=figure(2)
                imshow(ablationContour)
            hold on
            
                 
                horizontal = find(angleStore > horizontal_A -7 & angleStore <   horizontal_A +7);
               mid= round(length(horizontal)/2);
                
                     saveA_0 = boundary(horizontal(mid),2);
                     saveB_0 = boundary(horizontal(mid),1);
                     
                     plot(saveA_0,saveB_0,'xb','MarkerSize',12)
                     
              horizontal = find(angleStore >  horizontal_B - 7& angleStore < horizontal_B +7);
              midHorizontal = round(length(horizontal)/2);
              
                     saveA_1 = boundary(horizontal( midHorizontal ),2);
                     saveB_1 = boundary(horizontal( midHorizontal ),1);
                     
                     plot(saveA_1,saveB_1,'xb','MarkerSize',12);
                 
               vertical = find(angleStore >  horizontal_C -7 & angleStore <  horizontal_C +7);
               
              
               
                     saveA_2 = boundary(vertical(1),2);
                     saveB_2 = boundary(vertical(1),1);
                     
                     plot(saveA_2,saveB_2,'xr','MarkerSize',12);
                        verticalB = find(angleStore >  horizontal_D + -7 & angleStore <  horizontal_D +7);
                    
               verticalEnd = length(verticalB);
                     
                     saveA_3 = boundary(verticalB( verticalEnd),2);
                     saveB_3 = boundary(verticalB( verticalEnd),1);
                     
                     plot(saveA_3,saveB_3,'xr','MarkerSize',12);
                     
                       plot(STATS.Centroid(1),STATS.Centroid(2),'xr','MarkerSize',12);
                       
                        print(h,'-dpdf', ['Contour_EllipseIntersections',num2str(lauf),'.pdf']);%tifCurvetifFilename);
                     
                  if flag == 1
            
                   ecc_Store(lauf) = STATS.Eccentricity;
                   major_Store(lauf) = (pdist2([saveA_0(1) saveB_0(1)],[ saveA_1(1) saveB_1(1)])) * voxelX;
                   
                      minor_Store(lauf) =( pdist2([saveA_2(1) saveB_2(1)],[ saveA_3(1) saveB_3(1)])) * voxelX;
                   orientation_Store(lauf) = STATS.Orientation;
                   centroid_Store(lauf,1:2) = STATS.Centroid;
                    
                  else
                      
                    ecc_Store(lauf) = STATS.Eccentricity;
                    minor_Store(lauf) = (pdist2([saveA_0(1) saveB_0(1)],[ saveA_1(1) saveB_1(1)])) * voxelX;
                   major_Store(lauf)  =( pdist2([saveA_2(1) saveB_2(1)],[ saveA_3(1) saveB_3(1)])) * voxelX;
                    orientation_Store(lauf) = STATS.Orientation;
                    centroid_Store(lauf,1:2) = STATS.Centroid;
                      
                  end

            
            end
          %  figure(1)
             %  f
            tiffwrite_mat(imageFilterSave, [folderName,'_ROI'])  
            close all
          %  doDrawEllipse(STATS, STATS.Orientation,1,I_crop) 
            
            %clear boundary
            clear BW1
            if lauf < 10
                p = lauf -1
                
            else
                
            end
 
%%

for lauf =1:p
    clear xx yy
    
     cropOriginal =  imcrop(recoilImage(:,:,lauf),rect);
     imshow(cropOriginal,'InitialMagnification',1000)
     [xx,yy] = ginput(4)
     
       manualLongAxis = (pdist2([xx(1) yy(1)],[ xx(2) yy(2)])*voxelX)
       manualShortAxis = (pdist2([xx(3) yy(3)],[ xx(4) yy(4)])*voxelX)
       
       manualLongAxisStore(lauf) =  manualLongAxis;
       manualShortAxisStore(lauf) = manualShortAxis;
    
end
%p = lauf-1
%
%p=7
timeVec = ablationTime:timeInterval:(p+1)*timeInterval

%%%%% Identify chromation and calculate orientation of the cell

[cChromatin1 rChromatin1 cChromatin2 rChromatin2] = doChromatinChromatinDistance(imgMid)

         %%%% measurement 
          dx = cChromatin1(1) - cChromatin2(1);
          dy = rChromatin1(1) - rChromatin2(1);

          if dy < 0
              
               dx = cChromatin2(1) - cChromatin1(1);
               dy = rChromatin1(1) - rChromatin2(1);
                                   
                      angle = (atan2(dy,dx))*180/pi;
            
                  angle =    angle + 180        
                    

          else
               
               dx = cChromatin1(1) - cChromatin2(1);
               dy = rChromatin1(1) - rChromatin2(1);
                                   
              
                         angle = (atan2(dy,dx))*180/pi;
                         

                 
                     
          end
                     
                     
                     
                     
    
       
          
          %%%%%%% Assuming the the CR is orientated 90degree to the spindle
          %%%%%%% axis the 
          
          if angle < 0
              CR_angle = angle + 90
              
          else
            CR_angle = angle - 90
            
          end
          angle
         % CR_angle = angle +90

      
           %%%%%% Measure cleavage furrow vs pole ratio
           %%%%%% Mark first the polar regions and afterwards the
           %%%%%% contractile ring equator
           
        imshow(imgEquator(:,:,1),'InitialMagnification',1000)

           
           [xx,yy] = ginput(6)
           
           
           close all
           
           pole1 = (pdist2([xx(1) yy(1)],[ xx(2) yy(2)])) *midVoxelSize
           pole2 = (pdist2([xx(3) yy(3)],[ xx(4) yy(4)])) *midVoxelSize
           
           meanPole = (pole1+pole2)/2
           
           cleavageFurrow = (pdist2([xx(5) yy(5)],[ xx(6) yy(6)])) *midVoxelSize
           
           poleFurrowRatio = cleavageFurrow/meanPole
           
           
  %%%%%%%% Write data to disk     
  
    
            
            
            save([tifFilename,'.mat'])
            
          h=figure(1)   
            plot(timeVec,manualLongAxisStore(1:p),'ob')
            hold on
            plot(timeVec,major_Store,'or')
             hold on
       axis([0 10 0 10]) 
       set(gca,'FontSize',16,'FontName', 'Arial')

  %legend('Ra', 'SiR-actin pole')
  title('Outward movement after ablation - manual and automatic measurement long axis','FontSize', 16,'FontName', 'Arial')
    xlabel ('Time [s]','FontSize', 16,'FontName', 'Arial');
    ylabel('Outward movement after ablation [µm]' ,'FontSize', 16,'FontName', 'Arial');
       legend('Manual measurement', 'Automatic measurement')
  print(h,'-dpdf', [curdir '/' ,tifFilename,'_LifeactCherry_Recoil_automaticManual_LongAxis.pdf']);%tifCurvetifFilename);
  close all
  
  h=figure(1)   
            plot(timeVec,manualShortAxisStore(1:p),'ob')
            hold on
            plot(timeVec,minor_Store,'or')
             hold on
       axis([0 10 0 5]) 
       set(gca,'FontSize',16,'FontName', 'Arial')

  %legend('Ra', 'SiR-actin pole')
  title('Outward movement after ablation - manual and automatic measurement short axis','FontSize', 16,'FontName', 'Arial')
    xlabel ('Time [s]','FontSize', 16,'FontName', 'Arial');
    ylabel('Outward movement after ablation [µm]' ,'FontSize', 16,'FontName', 'Arial');
       legend('Manual measurement', 'Automatic measurement')
  print(h,'-dpdf', [curdir '/' ,tifFilename, '_LifeactCherry_Recoil_automaticManual_ShortAxis.pdf']);%tifCurvetifFilename);
   close all
  
 
   
   if  exist('orientation_Store') == 0

        cd([curdir]) 
  
            dummyVariable = (1:p)';
            dummyVariable(:) = 0;
           

            saveVariables = {};

            saveVariables{1} = timeVec(1:p)';%time{1};
            saveVariables{2} =   dummyVariable;
            
            saveVariables{3} =   dummyVariable;
            saveVariables{4} =   dummyVariable;
            saveVariables{5} =   dummyVariable;
            
          
            saveVariables{6} = dummyVariable;
           
           dummyVariable(1) = CR_angle;
            saveVariables{7} =  dummyVariable;
            
            dummyVariable(1) = meanPole;
            saveVariables{8} = dummyVariable;
            
            dummyVariable(1) = cleavageFurrow;
            saveVariables{9} = dummyVariable;
            
            
            dummyVariable(1) = poleFurrowRatio;
            saveVariables{10} = dummyVariable;
            
           
            
               dummyVariable(1) =midVoxelSize;
            saveVariables{11} = dummyVariable;
            
            
               dummyVariable(1) =voxelX;
            saveVariables{12} = dummyVariable;
            
             dummyVariable(1) =ablationTime;
            saveVariables{13} = dummyVariable;
            
             saveVariables{14} = manualLongAxisStore(1:p)';
             saveVariables{15} = manualShortAxisStore(1:p)';
            saveVariables{16} =      dummyVariable;
            saveVariables{17} =   dummyVariable;
            
            
           
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, ellipseOrientation, majorAxisEllipse, minorAxisEllipse',...
                ',eccentricityEllipse, orientationCell, orientationCleavageFurrow, meanDiameterAtPole',...
                ',meanDiameterFurrow, ratioFurrowPoleDiameter, midSectionVoxelSize , ablationVoxleSize, ablationTime',...
                ',manualMeasuredLongAxis , manualMeasuredShortAxis, contourRegularity, contourCircularity'];
            
            outid = fopen([tifFilename,'Analysis.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'Analysis.csv'],csvData,'roffset',1,'-append')
            
            
            %%%%% save only time and short axis
            
           

            saveVariables = {};

            saveVariables{1} = timeVec(1:p)';%time{1};
           
            
          
          
           
             saveVariables{2} = manualShortAxisStore(1:p)';
                saveVariables{3} =   dummyVariable;
            
           
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, shortAxisManual, shortAxisAutomatic' ];
            
            outid = fopen([tifFilename,'AnalysisSelected.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'AnalysisSelected.csv'],csvData,'roffset',1,'-append')
            
   else    
       
       
       
  cd([curdir]) 
  
            dummyVariable = (1:p)';
            dummyVariable(:) = 0;
           

            saveVariables = {};

            saveVariables{1} = timeVec(1:p)';%time{1};
            saveVariables{2} = orientation_Store(1:p)';
            
            saveVariables{3} = major_Store(1:p)';
            saveVariables{4} = minor_Store(1:p)';
            saveVariables{5} = ecc_Store(1:p)';
            
            dummyVariable(1) = angle;
            saveVariables{6} = dummyVariable;
           
           dummyVariable(1) = CR_angle;
            saveVariables{7} =  dummyVariable;
            
            dummyVariable(1) = meanPole;
            saveVariables{8} = dummyVariable;
            
            dummyVariable(1) = cleavageFurrow;
            saveVariables{9} = dummyVariable;
            
            
            dummyVariable(1) = poleFurrowRatio;
            saveVariables{10} = dummyVariable;
            
           
            
               dummyVariable(1) =midVoxelSize;
            saveVariables{11} = dummyVariable;
            
            
               dummyVariable(1) =voxelX;
            saveVariables{12} = dummyVariable;
            
             dummyVariable(1) =ablationTime;
            saveVariables{13} = dummyVariable;
            
             saveVariables{14} = manualLongAxisStore(1:p)';
             saveVariables{15} = manualShortAxisStore(1:p)';
            saveVariables{16} =     contourRegularity(1:p)';
            saveVariables{17} = contourCircularity(1:p)';
            
            
           
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, ellipseOrientation, majorAxisEllipse, minorAxisEllipse',...
                ',eccentricityEllipse, orientationCell, orientationCleavageFurrow, meanDiameterAtPole',...
                ',meanDiameterFurrow, ratioFurrowPoleDiameter, midSectionVoxelSize , ablationVoxleSize, ablationTime',...
                ',manualMeasuredLongAxis , manualMeasuredShortAxis, contourRegularity, contourCircularity'];
            
            outid = fopen([tifFilename,'Analysis.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'Analysis.csv'],csvData,'roffset',1,'-append')
            
            
            %%%%% save only time and short axis
            
           

            saveVariables = {};

            saveVariables{1} = timeVec(1:p)';%time{1};
           
            
          
          
           
             saveVariables{2} = manualShortAxisStore(1:p)';
                saveVariables{3} = major_Store(1:p)';
            
           
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, shortAxisManual, shortAxisAutomatic' ];
            
            outid = fopen([tifFilename,'AnalysisSelected.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'AnalysisSelected.csv'],csvData,'roffset',1,'-append')
            
            
   end


          

