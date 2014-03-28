function [Curvature] = doCalculateCurvatureAnisotropie(furrowRoi1_Sammel, slidingWindowSet, midFurrow1,greenStack,furrow1MidCoords)

slidingWindowSet = 40;

nn = length(furrowRoi1_Sammel);

for lauf = 1:nn
%%%% set sliding window for ellipse fitting


slidingWindow = slidingWindowSet;

%%%% crop the spline to a size as given by the sliding window

%InterSplineA_x = SplineXStore_1{lauf};
%InterSplineA_y = SplineYStore_1{lauf};


furrow1_Tmp  = round(furrowRoi1_Sammel{lauf})
%furrow2_Tmp  = round(furrowRoi2_Sammel{lauf})
 
    
runningSplineA_x = furrow1_Tmp(:,2);
runningSplineA_y = furrow1_Tmp(:,1);

%runningSplineB_x =  furrow2_Tmp(:,2);
%runningSplineB_y =  furrow2_Tmp(:,1);

InterSplineA_x = furrow1_Tmp(:,2);
InterSplineA_y = furrow1_Tmp(:,1);

%%%%%%%% Find the center pixel


centerA =  midFurrow1(lauf)

%find(rInterStore1(lauf) == InterSplineA_x & cInterStore1(lauf) == InterSplineA_y)


%centerA = find(cInterStore2(lauf) == test(:,2) & rInterStore2(lauf) == test(:,1))

if slidingWindow/2 > length(InterSplineA_x)-centerA

    slidingWindow = (length(InterSplineA_x)-centerA)/2

    else
end

splineA_x = InterSplineA_x(centerA-(slidingWindow)/2:centerA+(slidingWindow)/2)
splineA_y = InterSplineA_y(centerA-(slidingWindow)/2:centerA+(slidingWindow)/2)


InterSplineAStore{lauf} = [InterSplineA_y(centerA), InterSplineA_x(centerA)]
splineStoreA{lauf} = [splineA_x, splineA_y];

tmp = greenStack(:,:,1);
tmp(:,:) = 0;


%tmp(rNew1Ref{lauf}, cNew1Ref{lauf})=255;
%spline(rNew1Ref{lauf}', cNew1Ref{lauf}');

%plot(splineA_x,splineA_y,'-')

for sublauf = 1:length(InterSplineA_x)

    tmp(InterSplineA_x(sublauf),InterSplineA_y(sublauf)) = 255;
    
  %imshow(tmp,[])
  %  pause(0.1)
end

for sublauf = 1:length(splineA_x)
    tmp(splineA_x(sublauf),splineA_y(sublauf)) = 255;
end

    %imwrite(tmp,'curvatureTest','tif')
    
    circleParam = Kasa([splineA_y,splineA_x]);
    circleParam = LM([splineA_y,splineA_x],circleParam)
    
    imshow(tmp,[])
    hold on
    circle(circleParam(1),circleParam(2),circleParam(3))
    midCoordsTmp = furrow1MidCoords{lauf};
    plot(  midCoordsTmp(2), midCoordsTmp(1),'bx' )
    Curvature(lauf,:,:,:) = circleParam
close all

   

slidingWindow = slidingWindowSet;

end
   