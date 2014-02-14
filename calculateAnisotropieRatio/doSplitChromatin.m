function [B_FirstBoundary,B_SecondBoundary] = doSplitChromatin(ccc,BW2,flipMarker)

%%%% The function splits the image along the long axis of an ellipse fitted
%%%% around the segmented chromatin
k=1
tmp1 = doDrawLongAxis(ccc(k),BW2,flipMarker);

  

phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);
    
xbar = ccc(k).Centroid(1);
    ybar = ccc(k).Centroid(2);

    a = ccc(k).MajorAxisLength/2;
    b = ccc(k).MinorAxisLength/2;

    theta = pi*(ccc(k).Orientation/180);
    
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];

    xy = [a*cosphi; b*sinphi];
    xy = R*xy;

    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;

    imshow(BW2,[])
    hold on
    plot(x,y,'r','LineWidth',2);


%%%% To increase line thickness image could be filtered with a gaussian
    %G = fspecial('gaussian',[5 5],2);
  
    %imgFilter = imfilter(tmp1,G,'same');
    
    
    %%%%% Connect disconnected pixel
    tmp1 = bwmorph(tmp1,'bridge',8);
    
   
    %%%%%% identify individual objects within the image
    [B L N A] = bwboundaries(tmp1,8);
   
    %%%%%%% Split individual object. Object 1 is always the entire object,
    %%%%%%% the two second largest objects are chosen
     
    B_End =  length(B)
    
    %%%%%% Identify length of individual objects
    for subLauf = 1:B_End
        B_Tmp = B{subLauf};
        B_Sub(subLauf) = length(B_Tmp);
        
    end
    
    %%%% Sort identified objects, starting with the largest
    B_Sort = sort(B_Sub,'descend');
    
    %%%% Select second and third longest object
    B_FirstElement = B_Sort(2);
    B_SecondElement = B_Sort(3);
    
    B_FirstBoundaryIndex = find(B_Sub == B_FirstElement);
    B_SecondBoundaryIndex = find(B_Sub == B_SecondElement);
    
    
    %%%% If both objects have the same size, manual assignemnt of indices
    if length(B_FirstBoundaryIndex)> 1
        B_FirstBoundaryIndex = B_FirstBoundaryIndex(1);
        B_SecondBoundaryIndex = B_SecondBoundaryIndex(2);
    end
    
    %%%% If both objects have the same size, manual assignemnt of indices
    if length(B_SecondBoundaryIndex)>1
        B_FirstBoundaryIndex = B_FirstBoundaryIndex(1);
        B_SecondBoundaryIndex = B_SecondBoundaryIndex(2);
    end
    
    
    %%%% Store identified boundaries
    
    B_FirstBoundary = B{B_FirstBoundaryIndex};
    B_SecondBoundary = B{B_SecondBoundaryIndex};
    
   