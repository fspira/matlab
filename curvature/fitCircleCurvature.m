
 
 [mm pp nn] = size(greenStack);
tmp1 = greenStack(:,:,1);
nn = AnalysisEnd;
for lauf = 1:20%nn
tmp1(:,:)=0;
imshow(tmp1)
hold on
plot(cNew1Ref{lauf},rNew1Ref{lauf})
plot(cNew2Ref{lauf},rNew2Ref{lauf})
%plot(cInterStore1(lauf),rInterStore1(lauf),'xr')
%plot(cInterStore2(lauf),rInterStore2(lauf),'xr')
plot(cNewMid1(lauf),rNewMid1(lauf),'xy')
plot(cNewMid2(lauf),rNewMid2(lauf),'xy')
pause(0.1)
lauf
end
 

for lauf = 1:nn

   % midFurrow1(lauf,:) = [cInterStore1(lauf) rInterStore1(lauf)];
   % midFurrow2(lauf,:) = [cInterStore2(lauf) rInterStore2(lauf)];
     
     furrowRoi1_Sammel{lauf} = [cNew1Ref{:,lauf}  rNew1Ref{:,lauf}]
     furrowRoi2_Sammel{lauf} = [cNew2Ref{:,lauf}  rNew2Ref{:,lauf}]


end

%%%%%%%

%%%%%%%% circle fit test section


%%%%%% crop the roi for the circle fit
%%%%%%
%%%%%%
for lauf = 1:nn
slidingWindow = 35;
furrowRoi1TMP = furrowRoi1_Sammel{lauf};
furrowRoi2TMP = furrowRoi2_Sammel{lauf};



fitRegion1 = furrowRoi1TMP([redMax1(lauf)-slidingWindow:redMax1(lauf)+slidingWindow],1:2);

fitRegion2 = furrowRoi2TMP([redMax2(lauf)-slidingWindow:redMax2(lauf)+slidingWindow],1:2);


   greenRGB(:,:,:) = ind2rgb(normalizedImage(greenStack(:,:,lauf)),green);
   redRGB(:,:,:) = ind2rgb(normalizedImage(redStack(:,:,lauf)),red);
             
   tmpMerge(:,:,:,:) = greenRGB(:,:,:) + redRGB(:,:,:);

    %imwrite(tmp,'curvatureTest','tif')
    
  circleParam1 = Kasa([fitRegion1(:,1),fitRegion1(:,2)]);
 % circleParam1 = LM([fitRegion1(:,2),fitRegion1(:,1)], circleParam);
  
  
  circleParam2 = Kasa([fitRegion2(:,1),fitRegion2(:,2)]);
  %circleParam2 = LM([fitRegion2(:,2),fitRegion2(:,1)], circleParam);
    
    
    imshow(tmpMerge,[])
    hold on
    circle(circleParam1(1),circleParam1(2),circleParam1(3));
    circle(circleParam2(1),circleParam2(2),circleParam2(3));
 %   midCoordsTmp = furrow1MidCoords{lauf};
  %  plot(  centerA_Pos(1), centerA_Pos(2),'rx' )
    Curvature1(lauf,:,:,:) = circleParam1;
    Curvature1(lauf,:,:,:) = circleParam2;
%close all

region1_Store{lauf} = fitRegion1;
region2_Store{lauf} = fitRegion2;
   
%%%%%%%%
%%%% slidingWindowSet in µm
pause(0.2)

radiusSammelA(lauf) = circleParam1(3);
radiusSammelB(lauf) = circleParam2(3);
end


clear curvatureA CurvatureB CurvatureBNew CurvatureANew
clear curvatureB
clear CurvatureMean

CurvatureMaxA = find(max(radiusSammelA) == radiusSammelA);
CurvatureANew = cat(2,(0 - radiusSammelA(1:CurvatureMaxA)),radiusSammelA(CurvatureMaxA+1:length(radiusSammelA)));
 CurvatureA = (1./(CurvatureANew .* voxelX_mum) );
    
    
    
    CurvatureMaxB = find(max(radiusSammelB) == radiusSammelB);
     CurvatureBNew = cat(2,(0 - radiusSammelB(1:CurvatureMaxB)),radiusSammelB(CurvatureMaxB+1:end));
     CurvatureB = (1./(CurvatureBNew.* voxelX_mum));
    
     for lauf = 1:length(CurvatureA)
      
             CurvatureMean(lauf) =   (CurvatureA(lauf) +  CurvatureB(lauf))/2
    
     end
     
     CurvatureMean = CurvatureMean'
h=figure
     plot(timeVec(1:end), CurvatureMean,'x');
     axis([-200 +500 -0.4 1]) 
  %title('Gaussian fit to MRLIIC width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
   title('Curvature of the cleavage furrow. 2466, SiRActin MyoII','FontSize', 16)
    xlabel ('Time [s]','FontSize', 16);
    ylabel('Curvature [1/r]' ,'FontSize', 16);
     print(h,'-dpdf', [curdir , '/', tifFilename, '_Curvature.pdf']);%tifCurvetifFilename);
    
  %  legend([num2str(FWHMOrigStore(lauf)) 'µm' ' at ' num2str(FWHMTimeVec(lauf)) 's' ])

open CurvatureMean
  
            save('Workspace1FitCircle')
            
            timeVec = timeVec'

%close(gcf)

