function [Curvature,tmpMerge,InterSpline_Store] = doCalculateCurvatureAnisotropie(furrowRoi1_Sammel, slidingWindowSet, midFurrow1,midFurrow2,greenStack,green,red,centerA_Force)

%slidingWindowSet = 40;

%centerA_Pos =  midFurrow1(lauf,1:2); %%%% coordinates
%centerB_pos = midFurrow2(lauf,1:2);

nn = length(furrowRoi1_Sammel);
lauf = 1;
%for lauf = 1:nn
%%%% set sliding window for ellipse fitting


slidingWindow = slidingWindowSet;

%%%% crop the spline to a size as given by the sliding window

%InterSplineA_x = SplineXStore_1{lauf};
%InterSplineA_y = SplineYStore_1{lauf};


furrow1_Tmp  = round(furrowRoi1_Sammel{lauf});
%furrow2_Tmp  = round(furrowRoi2_Sammel{lauf})
 
    
runningSplineA_x = furrow1_Tmp(:,2);
runningSplineA_y = furrow1_Tmp(:,1);

%runningSplineB_x =  furrow2_Tmp(:,2);
%runningSplineB_y =  furrow2_Tmp(:,1);

InterSplineA_x = furrow1_Tmp(:,2);
InterSplineA_y = furrow1_Tmp(:,1);

%%%%%%%% Find the center pixel


centerA_Pos =  midFurrow1(1:2); %%%% coordinates
centerB_pos = midFurrow2(1:2);

%%%%%% find the position of the corresponding coordinates within the spline


 centerA = find(centerA_Pos(2) == InterSplineA_x & centerA_Pos(1) == InterSplineA_y);
 


 
 %%%% the furrows may sometimes be flipped, in this case use the center
 %%%% from the opposite spline
 
 if centerA == isempty(centerA)
      centerA = find(centerB_pos(2) == InterSplineA_x & centerB_pos(1) == InterSplineA_y);
 end
 
 [mCenter nCenter pCenter] = size(centerA);

    if mCenter > 1 

        if centerA(2) > 350
        centerA = centerA(1);
        else
            centerA  = centerA(2)
        end
    end
    
    %%%% if the correct midpoint can't be detected, its position can be enforced
    if centerA_Force == 0
    else
        centerA = centerA_Force;
        
    end
%centerA = find(cInterStore2(lauf) == test(:,2) & rInterStore2(lauf) == test(:,1))

if slidingWindow/2 > length(InterSplineA_x)-centerA

    slidingWindow = (length(InterSplineA_x)-centerA)/2;

    else
end

splineA_x = InterSplineA_x(centerA-(slidingWindow)/2:centerA+(slidingWindow)/2);
splineA_y = InterSplineA_y(centerA-(slidingWindow)/2:centerA+(slidingWindow)/2);


InterSplineAStore{lauf} = [InterSplineA_y(centerA), InterSplineA_x(centerA)];
splineStoreA{lauf} = [splineA_x, splineA_y];

tmp = greenStack(:,:,1);
tmp(:,:) = 0;

tmp1 = greenStack(:,:,1);
tmp1(:,:) = 0;

%tmp(rNew1Ref{lauf}, cNew1Ref{lauf})=255;
%spline(rNew1Ref{lauf}', cNew1Ref{lauf}');

%plot(splineA_x,splineA_y,'-')

for sublauf = 1:length(InterSplineA_x)

    tmp(InterSplineA_x(sublauf),InterSplineA_y(sublauf)) = 255;
    
  %imshow(tmp,[])
  %  pause(0.1)
end

for sublauf = 1:length(splineA_x)
    tmp1(splineA_x(sublauf),splineA_y(sublauf)) = 255;
end


   tmp1RGB(:,:,:) = ind2rgb(normalizedImage(tmp1(:,:)),green);
   tmpRGB(:,:,:) = ind2rgb(normalizedImage(tmp(:,:)),red);
             
   tmpMerge(:,:,:,:) = tmp1RGB(:,:,:) + tmpRGB(:,:,:);

    %imwrite(tmp,'curvatureTest','tif')
    
  circleParam = Kasa([splineA_y,splineA_x]);
    circleParam = LM([splineA_y,splineA_x], circleParam);
    
    imshow(tmpMerge,[])
    hold on
    circle(circleParam(1),circleParam(2),circleParam(3));
 %   midCoordsTmp = furrow1MidCoords{lauf};
    plot(  centerA_Pos(1), centerA_Pos(2),'rx' )
    Curvature(lauf,:,:,:) = circleParam;
%close all

InterSpline_Store = [splineA_x splineA_y];
   centerA

slidingWindow = slidingWindowSet;

%end
   