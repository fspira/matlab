function arrayOut = paddingWindow(arrayIn, slidingWindow)


if slidingWindow > 3
    
    
    placeHolder = zeros(1,((slidingWindow-1)/2));

    paddingArrayLeft = arrayIn(1) + placeHolder;
    paddingArrayRight = arrayIn(end) + placeHolder;
    
    paddingLeft = cat(2,paddingArrayLeft,arrayIn);
    paddingRight = cat(2,paddingLeft,paddingArrayRight);
    
    arrayOut = paddingRight;

elseif slidingWindow == 3
    
    
    placeHolder = 1;
    paddingArrayLeft = arrayIn(1) + placeHolder;
    paddingArrayRight = arrayIn(end) + placeHolder;
    
    paddingLeft = cat(2,paddingArrayLeft,arrayIn);
    paddingRight = cat(2,paddingLeft,paddingArrayRight);
    arrayOut = paddingRight;
    
elseif slidingWindow == 1
    
    arrayOut = arrayIn;
    
end


end