%%%%%%%% Function to detect curvature of the membrane

%[] = doDetermineCurvature()
[mm nn pp] = size(cNew1Ref);

greenStack = greenImg;
redStack = redImg;
 furrow1_Tmp  = furrowRoi1_Sammel{lauf}
 furrow2_Tmp  = furrowRoi2_Sammel{lauf}
 
 for lauf= 1:length(yFurrow2)
     rNewMid1(lauf) = round(yFurrow1{lauf});
     cNewMid1(lauf) = round(xFurrow1{lauf});
     
     rNewMid2(lauf) = round(yFurrow2{lauf});
     cNewMid2(lauf) = round(xFurrow2{lauf});

 end

for lauf = 1:nn
    close all
 furrow1_Tmp  = furrowRoi1_Sammel{lauf}
 furrow2_Tmp  = furrowRoi2_Sammel{lauf}
 
    
runningSplineA_x = furrow1_Tmp(:,1);
runningSplineA_y = furrow1_Tmp(:,2);

runningSplineB_x =  furrow2_Tmp(:,1);
runningSplineB_y =  furrow2_Tmp(:,2);

%%%% identifies the central region of the spline
%center =   redMax1(lauf)
tmp1 = greenStack(:,:,1);
tmp1(:,:) = 0;

tmp1(rNewMid1(lauf), cNewMid1(lauf)) =255;

tmp = greenStack(:,:,1);
tmp(:,:) = 0;
imshow(tmp1)
hold on
line(runningSplineB_x,runningSplineB_y)

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

yCat = [y1' [1:512]']

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



for subrun = 1:length(runningSplineA_x)
    tmpSpline(round(runningSplineA_y(subrun)),round(runningSplineA_x(subrun))) = 125;

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



for subrun = 1:length(runningSplineB_x)
    tmpSpline(round(runningSplineB_y(subrun)),round(runningSplineB_x(subrun))) = 125;

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
pause(0.1)
lauf
end
 

[Curvature] = doCalculateCurvatureAnisotropie(furrowRoi1_Sammel, slidingWindowSet, midFurrow1,greenStack)




L2 = greenStack(:,:,1);
L2(:,:) = 0;
for lauf = 1:nn
    
    
  imshow(L2,[])
    hold on
 %line(SplineYStore_1{lauf},  SplineXStore_1{lauf})
 splineStoreTmp = splineStoreA{lauf}
 line(splineStoreTmp(:,2),splineStoreTmp(:,1))
  % plot(flank1Store{:lauf},flank1Store{lauf}, 'b', 'LineWidth', 2)
  
  circleTmp =Curvature{lauf}
  
  %circle(circleTmp(2),circleTmp(1),circleTmp(3))
  
 centerTmp =  InterSplineAStore{lauf}
  
  %  plot( cNew1Ref{lauf}, rNew1Ref{lauf}, 'g', 'LineWidth', 2)
    plot(centerTmp(:,1),centerTmp(:,2) ,'xr')
    
    pause(0.2)
     %title(['MMLCII-B Furrow Region:' tifFilename],'FontSize', 20);
   % mov(lauf) = getframe(gcf);
     hold off
     
     radiusSammel(lauf) = circleTmp(3)
end

plot(1:lastFrameToConsider,radiusSammel,'x')
hold on

x1 = [1:lastFrameToConsider]';
x2 = radiusSammel'
k=10

[p,S,mu] = polyfit(x1, x2, k);
    
    %%%%% Test the model by evaluating test data. The test data is real
    %%%%% measured data
    %%%%% Predicted time line is the fitted curve. Data points for staging
    %%%%% will be compared against this time line.
    predicted_Curvature = polyval(p, x2); 
    
    plot(1:lastFrameToConsider,predicted_Curvature,'-')


close(gcf)

%# save as AVI file, and open it using system video player
movie2avi(mov, 'Ellipse.avi', 'compression','None', 'fps',10);
 
cNewMid1(lauf), rNewMid1(lauf)

spline(cNew1Ref{lauf}', rNew1Ref{lauf}')
