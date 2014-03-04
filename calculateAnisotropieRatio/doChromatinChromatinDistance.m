function [cChromatin1 rChromatin1 cChromatin2 rChromatin2] = doChromatinChromatinDistance(imgMid)

            
   %MIJ.setRoi( [ cChromatin1(lauf); rChromatin1(lauf)], ij.gui.Roi.POINT);
   %             k = waitforbuttonpress 
   %             coords =    MIJ.getRoi(1);
                   
    %          cChromatin1(lauf) = coords(1);
    %          rChromatin1(lauf) = coords(2);
            
    %              MIJ.setRoi( [ cChromatin2(lauf); rChromatin2(lauf)], ij.gui.Roi.POINT);
    %            k = waitforbuttonpress 
    %            coords =    MIJ.getRoi(1);
                   
    %           cChromatin2(lauf) = coords(1);
    %            rChromatin2(lauf) = coords(2);

    %            MIJ.run('closeAllWindows');
%imgMid = greenStack;
[mm nn pp] = size(imgMid)
AxisDiam = 1;


  flipMarker = 0;
  ellipseOrientation = [];
  ellipseOrientationStore =[];
  chromatin1 ={};
  chromatin2 = {};
  
for lauf= 1:pp
    
   % imshow(imgMid(:,:,lauf),[])
    %hold on
    %plot(center(1),center(2),'x')
   
     
    G1 = fspecial('gaussian',[75 75],25);
    imgFilter = imfilter(imgMid(:,:,lauf),G1,'same');
    BW2 = edge(imgFilter,'canny',0.7,10);
    BW2 = bwmorph(BW2,'bridge',8);
    imshow(BW2)
    
    ccc = regionprops(BW2,'all');
    [structLength, tmp] = size(ccc)
    axisRatio = ccc(1).MajorAxisLength/ccc(1).MinorAxisLength
  
    imshow(BW2)
    
    
    %%%%% This section was implemented to correctly segment early
    %%%%% anaphases. Chromatin may be connected but 95% percent of the mass
    %%%%% is separated. A Big smoothing will identify a single object,
    %%%%% which can be correctly split.
    
    if axisRatio < 1.2 && structLength > 1
            
            imshow(imgFilter,[])
            G = fspecial('gaussian',[80 80],25);
            imgFilter = imfilter(imgMid(:,:,lauf),G,'same');
            BW2 = edge(imgFilter,'canny',0.4,10);
            ccc = regionprops(BW2,'all');
            structLength = length(ccc);
    end
    
    %%%%% This section splits a single object into two. A large gaussian is
    %%%%% applied to smooth out the objects. Decision is made upon the axis
    %%%%% ratio (not elongated any more and the number of identified
    %%%%% objects.
    
    if axisRatio < 1.6 && structLength == 1
  
            imshow(imgFilter,[])
            G = fspecial('gaussian',[80 80],25);
            imgFilter = imfilter(imgMid(:,:,lauf),G,'same');
            BW2 = edge(imgFilter,'canny',0.4,10);
            ccc = regionprops(BW2,'all');


            %k=1
            imshow(BW2)
            %hold on
            %plot(ccc(k).Centroid(1),ccc(k).Centroid(2),'xr')
            %imshow(imgMid(:,:,lauf),[])


            %%%%%%% Dilation and erosion is required to close gaps in the outline
          %  
           %se = strel('disk',6);

          % BW3 = imdilate(BW2,se);

           % BW4 = imerode(BW3,se,6);
           % BW4 = BW2;


            %%%%% Deterimine the ellipse orientation 1 = vertical 0 =  horizontal

                    if abs(ccc(1).Orientation) > 40 


                        ellipseOrientation = 1;


                    else


                        ellipseOrientation = 0;


                    end

                    %%%%%% Test whether ellipse has been flipped, if ellipse is
                    %%%%%% flipped between consecutive frames the old orientation
                    %%%%%% will be used. This anticipates that chromatin orientation should
                    %%%%%% not change during cytokinesis.

                    varCheck = isempty(ellipseOrientationStore);

                    if varCheck == 0 
                        if ellipseOrientation == ellipseOrientationStore

                        else

                            ellipseOrientation = ellipseOrientationStore
                            ccc(1).Orientation = orientationStore
                            flipMarker = 1
                        end
                    end




            %plot(ccc(1).Centroid(1),ccc(1).Centroid(2),'x')
            %plot(ccc(2).Centroid(1),ccc(2).Centroid(2),'x')

           %  if axisRatio < 1.6

                    %%%%% split chromatin along the long axis of the ellipse



                    [B_FirstBoundary,B_SecondBoundary] = doSplitChromatin(ccc(1),BW2,flipMarker);
                    [first_css, second_css] = doCalculateCenter(B_FirstBoundary,B_SecondBoundary,BW2,imgMid(:,:,lauf))
                      close all

                    h= figure(1)
                    imshow(imgMid(:,:,lauf),[])
                    hold on
                    plot(first_css.WeightedCentroid(1),first_css.WeightedCentroid(2),'xr')
                    plot(second_css.WeightedCentroid(1),second_css.WeightedCentroid(2),'xb') 
            
            chromatin1{lauf} = first_css.WeightedCentroid;
         
            chromatin2{lauf} = second_css.WeightedCentroid;
            
            


                          print(h,'-dpdf', [num2str(lauf),'_Chromatin.pdf']);  

                    close all

                      orientationStore = ccc.Orientation;
                      centerStore = ccc.Centroid;


                         orientationStore = ccc(1).Orientation;
                         ellipseOrientationStore = ellipseOrientation;

                         flipMarker = 0
                         
    else
        %%%%% Orientation of the chromatin is identified and stored.
        %%%%% Orientations are grouped to horizontal and vertical.
        
                    if abs(ccc(1).Orientation) > 40 


                        ellipseOrientation = 1;


                    else


                        ellipseOrientation = 0;


                    end

                    %%%%%% Test whether ellipse has been flipped, if ellipse is
                    %%%%%% flipped between consecutive frames the old orientation
                    %%%%%% will be used. This anticipates that chromatin orientation should
                    %%%%%% not change during cytokinesis.

                    varCheck = isempty(ellipseOrientationStore);

                    if varCheck == 0 
                        if ellipseOrientation == ellipseOrientationStore

                        else

                            ellipseOrientation = ellipseOrientationStore
                            ccc(1).Orientation = orientationStore
                            flipMarker = 1
                        end
                    end




            %plot(ccc(1).Centroid(1),ccc(1).Centroid(2),'x')
            %plot(ccc(2).Centroid(1),ccc(2).Centroid(2),'x')

           %  if axisRatio < 1.6

                    %%%%% split chromatin along the long axis of the ellipse

                  BW2 = bwmorph(BW2,'bridge',8);
                  
                  
                  s = regionprops(BW2, imgMid(:,:,lauf), {'WeightedCentroid'});
   
if structLength == 1
                    h= figure(1)
                    imshow(imgMid(:,:,lauf),[])
                    hold on
                    plot(s(1).WeightedCentroid(1),s(1).WeightedCentroid(2),'xb')
                   
                       %pause(0.5) 


                          print(h,'-dpdf', [num2str(lauf),'_Chromatin.pdf']);  

                    close all
                              
            chromatin1{lauf} = s(1).WeightedCentroid;
         
            chromatin2{lauf} = s(1).WeightedCentroid;
else
    
             h= figure(1)
                    imshow(imgMid(:,:,lauf),[])
                    hold on
                    plot(s(1).WeightedCentroid(1),s(1).WeightedCentroid(2),'xr')
                    plot(s(2).WeightedCentroid(1),s(2).WeightedCentroid(2),'xb')
                   
                       %pause(0.5) 


                          print(h,'-dpdf', [num2str(lauf),'_Chromatin.pdf']);  

                    close all
                    
                    
            chromatin1{lauf} = s(1).WeightedCentroid;
         
            chromatin2{lauf} = s(2).WeightedCentroid;
                    
end

                      orientationStore = ccc.Orientation;
                      centerStore = ccc.Centroid;


                         orientationStore = ccc(1).Orientation;
                         ellipseOrientationStore = ellipseOrientation;

                         flipMarker = 0
                         

    
    end
end

   
    

    
    %%%%% Compute the weighted center of mass and connect uncontected
    %%%%% pixels before computing the center of mass
  
   
   
   
  %  imshow(BW4)
  %  hold on
  %  plot(s(1).Centroid(1),s(1).Centroid(2),'xr')
  %  plot(s(2).Centroid(1),s(2).Centroid(2),'xr')
  %  plot(chromatin1(1),chromatin1(2),'xb')
  %   plot(chromatin2(1),chromatin2(2),'xb')
  %  line([chromatin1(1) chromatin2(1)],[chromatin1(2) chromatin2(2)])
     
   % close all
    
           for subrun = 1:length(chromatin1)

           % Calculate the angle between the points and convert it from rad to the
           % degree system

             % dx = chromatin2(1) - chromatin1(1);
             % dy = chromatin2(2) - chromatin1(2);

             % angle = (abs(atan2(dy,dx)))*180/pi




              % BW_test = BW2;
             %  BW_test(:,:) = 0;
          %    imshow(BW_test)
                %hold on


                % dist2 = pdist2([ chromatin1(1),chromatin1(2)],[ chromatin2(1), chromatin2(2)]);

                 %distChrom = dist2*voxelX_mumMid;

                chromatinTmp = chromatin1{subrun}
                cChromatin1(subrun) = chromatinTmp(:,1);
                
                rChromatin1(subrun) = chromatinTmp(:,2);

                chromatinTmp = chromatin2{subrun}
                cChromatin2(subrun) = chromatinTmp(:,1);
                
                rChromatin2(subrun) = chromatinTmp(:,2);
           end
    end

    