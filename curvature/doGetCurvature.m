
function [Curvature,tmp]= doGetCurvature(rIntersect2,cIntersect2,splineFitOutline,slidingWindow,BW1)
lauf = 1;

%slidingWindow = 30;


centerA = find(rIntersect2(1) == splineFitOutline(:,1) & cIntersect2(1) ==  splineFitOutline(:,2));

InterSplineA_x = splineFitOutline(:,1);
InterSplineA_y = splineFitOutline(:,2);


splineA_x = InterSplineA_x(centerA-(slidingWindow)/2:centerA+(slidingWindow)/2);

splineA_y = InterSplineA_y(centerA-(slidingWindow)/2:centerA+(slidingWindow)/2);


InterSplineAStore{lauf} = [InterSplineA_y(centerA), InterSplineA_x(centerA)];
splineStoreA{lauf} = [splineA_x, splineA_y];

tmp = BW1(:,:,1);
tmp(:,:) = 0;


%tmp(rNew1Ref{lauf}, cNew1Ref{lauf})=255;
%spline(rNew1Ref{lauf}', cNew1Ref{lauf}');

%plot(splineA_x,splineA_y,'-')

for sublauf = 1:length(splineA_x)

    tmp(splineA_x(sublauf),splineA_y(sublauf)) = 255;
   % imshow(tmp,[])
   % pause(0.1)
end


    %imwrite(tmp,'curvatureTest','tif')
    
    circleParam = Kasa([splineA_y,splineA_x]);
    circleParam = LM([splineA_y,splineA_x],circleParam);
    
    h= figure
    imshow(tmp,[])
    hold on
    circle(circleParam(1),circleParam(2),circleParam(3));
    
    
    Curvature = circleParam;