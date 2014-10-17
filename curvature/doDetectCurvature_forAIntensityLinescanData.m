%%%%%%%% Function to detect curvature of the membrane

%[] = doDetermineCurvature()



distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;

blue = zeros(256,3);
blue(:,3) = distvec;

curdir = pwd;


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


%Miji;

%[m n p] = size(greenImg);




  
  [mm nn pp] = size(cNew1Ref);
 
  for lauf = 1:nn
  
    flank1NewStore{lauf} =  [ cNew1Ref{lauf} rNew1Ref{lauf}]
    flank2NewStore{lauf} =  [ cNew2Ref{lauf} rNew2Ref{lauf}]

  end
 
for lauf = 1:nn
    close all
 furrow1_Tmp  =  flank1NewStore{lauf} %flank1Store{lauf}
 furrow2_Tmp  = flank2NewStore{lauf}%flank2Store{lauf}
 

 className = class(furrow1_Tmp);
 
 if strcmp(className ,'double')
    runningSplineA_x = furrow1_Tmp(:,1);
    runningSplineA_y = furrow1_Tmp(:,2);

    runningSplineB_x =  furrow2_Tmp(:,1);
    runningSplineB_y =  furrow2_Tmp(:,2);
    
 else

    runningSplineA_x = furrow1_Tmp{:,1};
    runningSplineA_y = furrow1_Tmp{:,2};

    runningSplineB_x =  furrow2_Tmp{:,1};
    runningSplineB_y =  furrow2_Tmp{:,2};

 end
%%%% identifies the central region of the spline
%center =   redMax1(lauf)
tmp1 = greenStack(:,:,1);
tmp1(:,:) = 0;

tmp1(rNewMid1(lauf), cNewMid1(lauf)) =255;

tmp = greenStack(:,:,1);
tmp(:,:) = 0;
imshow(tmp1)
hold on
line(runningSplineA_x,runningSplineA_y)

plot(cNewMid1(lauf),rNewMid1(lauf),'xb')
plot(cNewMid2(lauf),rNewMid2(lauf),'xr')
pause(0.1)



%%%% Calculate intersection points between ingressing furrow
%%%% First exrapolate line between marked center points of the ingressing
%%%% furrow
imshow(tmp)
xp = [cNewMid1(lauf) cNewMid2(lauf)];
yp = [rNewMid1(lauf) rNewMid2(lauf)];


    a = (yp(2)-yp(1)) / (xp(2)-xp(1));
    b = yp(1)-a*xp(1);

  %%%% Determine size of the active window
  xlims = xlim(gca);
    ylims = ylim(gca);

%%%% Calculate line
y = xlims*a+b;

hold on
plot(xp(1),yp(1),'xb')
plot(xp(2),yp(2),'xr')
line( xlims, y );

for subrun = 1:512
    
    y1(subrun) = subrun*a+b;

end

y1shape = y1(round(xlims(1))):round(xlims(2))
    
%%%% store values of the line connecting the center points of the
%%%% ingressing furrow
%%%% These lines eliminate negative values
%yCat = [y1' [round(xlims(1)):round(xlims(2))]']

yCat = [y1 [1:512]']

yCatTmp = yCat(:,1:2) < xlims(2)
yCat = yCat .* yCatTmp;
yCat = [round(yCat(:,1)),round(yCat(:,2))];

yCatNorm = (yCat > 0);
yCat = yCatNorm .* yCat;
catZero = find(yCat(:,1) >0);

yCat = yCat(catZero,:);
%catNew = yCat(catZero)

%flag =  isempty(catZero)
%if flag == 0 

 %   yCat= yCat(1:catZero(1)-1,1:2);
 %   else
%end

%%%% To obtain values for intersecting points geometrically draw the shape
%%%% and the line connecting two points assigning values of 125 to each
%%%% shape. After addition of both shapes pixel values having a value
%%%% greate than 125 are intersecting pixels.
%%%% Calculate flank 1
tmp(:,:) = 0;
tmp = double(tmp);
tmp1 = tmp;
tmp1(:,:) = 0;
knots = [];

%%%%% This section fits a spline to the individual data points
 splineEnd = length(runningSplineA_x);
 
 knots(:,1:splineEnd) = [runningSplineA_x runningSplineA_y]';
%x = knots(1, :);
%y = knots(2, :);
%areaOfPolygon = polyarea(x,y);
finerSpacing = 1:0.01:splineEnd;
originalSpacing = 1:length(knots);

splineXY = spline(originalSpacing, knots, finerSpacing);

A = round(splineXY);
splineXY = A';

%%%% The following tow lines will reduce the number of knots within the
%%%% spline. 

%%[u,id1,id2] = unique(A(:,1:2),'rows');

%%splineXY = u';


%plot(knots(1, :), knots(2, :), 'ro',...
%splineXY(1, :), splineXY(2, :), 'b+-', 'LineWidth', 2, 'MarkerSize',16);

imshow(tmp)
hold on
tmpSpline = tmp;

X1 = splineXY(:,1);
Y1 = splineXY(:,2);



for subrun = 1:length(X1)
    tmpSpline(round(Y1(subrun)),round(X1(subrun))) = 125;

end
imshow(tmpSpline,[])

tmp= tmpSpline;

%%%%%%%%% This section re-initializes the index of the BW image. This is
%%%%%%%%% important because the spline function uses a left to right
%%%%%%%%% indexing


B = bwboundaries(tmp,8);

 BLength = length(B);


[m n p] = size(B)

for subrun = 1:m
    
    BLength(subrun) = length(B{subrun});
    
end

B_Select = max(BLength);
B_Select = find(B_Select == BLength);
B = B{B_Select};
SplineXStore_1{lauf} = B(:,1);
SplineYStore_1{lauf} =  B(:,2);

tmp(:,:) = 0;

for subrun  = 1:length(B(:,1))
    tmp(B(subrun,1),B(subrun,2)) = 125;
end

imshow(tmp,[])

clear B
     
%%%%% this section plots xy coordinates marking ingressing furrow into the
%%%%% image and determines the intersection spline coordinates

     close all
     imshow(tmp)
     hold on
     plot(xp(1),yp(1),'xr')
     plot( xp(2),yp(2),'xr')
     hLine = imline(gca,[xp(1),xp(2)], [yp(1),yp(2)]); 
     tmp1 = hLine.createMask();
     close all

     
tmp1 = im2bw(tmp1);
tmp1 = tmp1.*125;

tmp2 = tmp + tmp1;
imshow(tmp2)

[rIntersect1,cIntersect1] = find(tmp2 > 130)

flag = isempty(rIntersect1)

        
        if flag == 1
                    [yCat] = doInterpolateLineThroughPoints(xp, yp, tmp)
                       
              test = yCat(yCat(:,2)~=0);
              yCat = yCat(1:length(test),1:2)

                for subrun = 1:length(yCat)

                    tmp1((round(yCat(subrun,1))),round(yCat(subrun,2))) = 125;

                end



                tmp1 = bwmorph(tmp1,'bridge',8);
               
                tmp1 = im2bw(tmp1);
                tmp1 = tmp1.*125;
                tmp1 = doAGaussianFiltering(tmp1);

                tmp2 = tmp + tmp1;



                [rIntersect1,cIntersect1] = find(tmp2 > 125);
                
          
            rIntersect1 = rIntersect1(1);
            cIntersect1 = cIntersect1(1);

                else
        end



%%%% use central value of line
[mi ni pi] = size(rIntersect1);
if mi >1
    rIntersect1 = rIntersect1(round(mi/2));
    cIntersect1 = cIntersect1(round(mi/2));
    else
end
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% calculate flank 2 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp(:,:) = 0;


%%%%% This section fits a spline to the individual data points
 splineEnd = length(runningSplineB_x);
 
 knots(:,1:splineEnd) = [runningSplineB_x runningSplineB_y]';
%x = knots(1, :);
%y = knots(2, :);
%areaOfPolygon = polyarea(x,y);
finerSpacing = 1:0.01:splineEnd;
originalSpacing = 1:length(knots);


splineXY = spline(originalSpacing, knots, finerSpacing);
%plot(knots(1, :), knots(2, :), 'ro',...
%splineXY(1, :), splineXY(2, :), 'b+-', 'LineWidth', 2, 'MarkerSize',16);


A = round(splineXY)';
splineXY= A;
%%%% The following tow lines will reduce the number of knots within the
%%%% spline. 
%[u,id1,id2] = unique(A(:,1:2),'rows');

%splineXY = u';


imshow(tmp,[])
hold on
tmpSpline = tmp;

imshow(tmpSpline,[])

X1 = splineXY(:,1);
Y1 = splineXY(:,2);



for subrun = 1:length(X1)
    tmpSpline(round(Y1(subrun)),round(X1(subrun))) = 125;

end
imshow(tmpSpline)

tmp= tmpSpline;

%%%%%%%%% This section re-initializes the index of the BW image. This is
%%%%%%%%% important because the spline function uses a left to right
%%%%%%%%% indexing

B = bwboundaries(tmp,8);


 BLength = length(B);

[m n p] = size(B)
for subrun = 1:m
    
    BLength(subrun) = length(B{subrun});
    
end

B_Select = max(BLength);
B_Select = find(B_Select == BLength);
B = B{B_Select};

SplineXStore_2{lauf} =  B(:,1);
SplineYStore_2{lauf} =  B(:,2);

tmp(:,:) = 0;

for subrun  = 1:length(B(:,1))
    tmp(B(subrun,1),B(subrun,2)) = 125;
end

clear B
%se = strel('ball',3,3,8); % radius of 4
%BW1 = imclose(tmp1,se);



imshow(tmp,[])
      
     close all
     imshow(tmp)
     hold on
     plot(xp(1),yp(1),'xr')
     plot( xp(2),yp(2),'xr')
     hLine = imline(gca,[xp(1),xp(2)], [yp(1),yp(2)]); 
     tmp1 = hLine.createMask();
     
     %tmp1 = uint8(tmp1) + uint8( binaryImage1);

     
    tmp1 = im2bw(tmp1);
    tmp1 = tmp1.*125;

    tmp2 = tmp + tmp1;
    imshow(tmp2)


[rIntersect2,cIntersect2] = find(tmp2 > 135)

flag = isempty(rIntersect2)

    
    if flag == 1
                [yCat] = doInterpolateLineThroughPoints(xp, yp, tmp)

                
              test = yCat(yCat(:,2)~=0);
              yCat = yCat(1:length(test),1:2)
                
            for subrun = 1:length(yCat)

                tmp1((round(yCat(subrun,1))),round(yCat(subrun,2))) = 125;

            end


        
            tmp1 = bwmorph(tmp1,'bridge',8);
            tmp1 = im2bw(tmp1);
            tmp1 = tmp1.*125;
            tmp1 = doAGaussianFiltering(tmp1);
           

            tmp2 = tmp + tmp1;


            [rIntersect2,cIntersect2] = find(tmp2 > 135);
            rIntersect2 = rIntersect2(1);
            cIntersect2 = cIntersect2(1);

            else
    end

[mi ni pi] = size(rIntersect2);

if mi > 1
    rIntersect2 = rIntersect2(round(mi/2));
    cIntersect2 = cIntersect2(round(mi/2));
    else
end
close all

cInterStore1(lauf) = cIntersect1(1);
rInterStore1(lauf) = rIntersect1(1);

cInterStore2(lauf) = cIntersect2(1);
rInterStore2(lauf) = rIntersect2(1);

close all

end




for lauf = 1:nn
tmp1(:,:)=0;
imshow(tmp1)
hold on
plot(cNew1Ref{lauf},rNew1Ref{lauf})
plot(cNew2Ref{lauf},rNew2Ref{lauf})
plot(cInterStore1(lauf),rInterStore1(lauf),'xr')
plot(cInterStore2(lauf),rInterStore2(lauf),'xr')
plot(cNewMid1(lauf),rNewMid1(lauf),'xy')
plot(cNewMid2(lauf),rNewMid2(lauf),'xy')
%pause(0.1)
lauf
end
 

for lauf = 1:nn

     midFurrow1(lauf,:) = [cInterStore1(lauf) rInterStore1(lauf)];
     midFurrow2(lauf,:) = [cInterStore2(lauf) rInterStore2(lauf)];
     
     furrowRoi1_Sammel{lauf} = [SplineYStore_1{:,lauf}  SplineXStore_1{:,lauf}]
     furrowRoi2_Sammel{lauf} = [SplineYStore_2{:,lauf}  SplineXStore_2{:,lauf}]


end


%%%% slidingWindowSet in µm
slidingWindowSet =3;


%lauf = 2

for lauf = 1:nn
   % slidingWindowSet = 40;


    [Curvature,tmpMerge,InterSpline_Store] = doCalculateCurvatureAnisotropie(furrowRoi1_Sammel(lauf), slidingWindowSet, midFurrow1(lauf,1:2),midFurrow2(lauf,1:2),greenStack(:,:,lauf),green,red,0,voxelX_mum);

    
    
    tmpMergeSammel_A(:,:,:,lauf) = tmpMerge;
    CurvatureSammel_A{lauf} = Curvature;
 
    [Curvature,tmpMerge,InterSpline_Store] = doCalculateCurvatureAnisotropie(furrowRoi2_Sammel(lauf), slidingWindowSet, midFurrow2(lauf,1:2),midFurrow1(lauf,1:2),greenStack(:,:,lauf),green,red,0,voxelX_mum);

    
    tmpMergeSammel_B(:,:,:,lauf) = tmpMerge;
    CurvatureSammel_B{lauf} = Curvature;
    InterSpline_Store_Sammel{lauf} = InterSpline_Store;

end

tmpMergeMerge = tmpMergeSammel_A + tmpMergeSammel_B;

for lauf = 1:nn
   
    %imshow(imgMerge(:,:,:,lauf),[])
    imshow(tmpMergeMerge(:,:,:,lauf),[])
    
   % pause(0.1)
    
    
end

L2 = greenStack(:,:,1);
L2(:,:) = 0;
for lauf = 1:nn
    
    
  %imshow(L2,[])
   imshow(tmpMergeMerge(:,:,:,lauf),[])
    hold on
 %line(SplineYStore_1{lauf},  SplineXStore_1{lauf})
 splineStoreTmpA =  furrowRoi1_Sammel{:,lauf};
 splineStoreTmpB =  furrowRoi2_Sammel{:,lauf};
 
 %line(splineStoreTmpA(:,1),splineStoreTmpA(:,2),'LineWidth',2,'Color','r')
 %line(splineStoreTmpB(:,1),splineStoreTmpB(:,2),'LineWidth',2,'Color','r')
 % plot(flank1Store{:lauf},flank1Store{lauf}, 'b', 'LineWidth', 2)
  
  circleTmpA =CurvatureSammel_A{lauf}
  circleTmpB =CurvatureSammel_B{lauf}
  
 circle(circleTmpA(1),circleTmpA(2),circleTmpA(3))
  
 circle(circleTmpB(1),circleTmpB(2),circleTmpB(3))
  
  
  
% centerTmp =  InterSplineAStore{lauf}
  
   %plot( cNew1Ref{lauf}, rNew1Ref{lauf}, 'g', 'LineWidth', 2)
 %  plot(centerTmp(:,2),centerTmp(:,1) ,'xr')
    
    pause(0.2)
     %title(['MMLCII-B Furrow Region:' tifFilename],'FontSize', 20);
    mov(lauf) = getframe(gcf);
     hold off
     
     radiusSammelA(lauf) = circleTmpA(3);
     radiusSammelB(lauf) = circleTmpB(3);
     
    
    
end

clear curvatureA CurvatureB CurvatureBNew CurvatureANew
clear curvatureB
clear CurvatureMean

CurvatureMaxA = find(max(radiusSammelA) == radiusSammelA);
CurvatureANew = cat(2,(0 - radiusSammelA(1:CurvatureMaxA)),radiusSammelA(CurvatureMaxA+1:length(radiusSammelA)));
 CurvatureA = (1./CurvatureANew);
    
    
    
    CurvatureMaxB = find(max(radiusSammelB) == radiusSammelB);
     CurvatureBNew = cat(2,(0 - radiusSammelB(1:CurvatureMaxB)),radiusSammelB(CurvatureMaxB+1:end));
     CurvatureB = (1./CurvatureBNew);
    
     for lauf = 1:length(CurvatureA)
      
             CurvatureMean(lauf) =   (CurvatureA(lauf) +  CurvatureB(lauf))/2
    
     end
h=figure
     plot(timeVec(1:end), CurvatureMean,'x');
     axis([-200 +500 -0.05 0.1]) 
  %title('Gaussian fit to MRLIIC width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
   title('Curvature of the cleavage furrow. 2466, SiRActin MyoII','FontSize', 16)
    xlabel ('Time [s]','FontSize', 16);
    ylabel('Curvature [1/r]' ,'FontSize', 16);
     print(h,'-dpdf', [curdir , '/', tifFilename, '_Curvature.pdf']);%tifCurvetifFilename);
    
  %  legend([num2str(FWHMOrigStore(lauf)) 'µm' ' at ' num2str(FWHMTimeVec(lauf)) 's' ])

CurvatureMean = CurvatureMean';

%close(gcf)

%movie2avi(mov, 'AnalysisROI.avi', 'compression','None', 'fps',10,'quality',100);


%plot(1:lastFrameToConsider,radiusSammel,'x')
%hold on

%x1 = [1:lastFrameToConsider]';
%x2 = radiusSammel'
%k=10

%[p,S,mu] = polyfit(x1, x2, k);
    
    %%%%% Test the model by evaluating test data. The test data is real
    %%%%% measured data
    %%%%% Predicted time line is the fitted curve. Data points for staging
    %%%%% will be compared against this time line.
   % predicted_Curvature = polyval(p, x2); 
  %  
   % plot(1:lastFrameToConsider,predicted_Curvature,'-')




%# save as AVI file, and open it using system video player
 
%cNewMid1(lauf), rNewMid1(lauf)

%spline(cNew1Ref{lauf}', rNew1Ref{lauf}')
