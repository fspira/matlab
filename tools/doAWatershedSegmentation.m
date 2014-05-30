
function[t3Store] = doAWatershedSegmentation(greenStack, redStack)


%%%%%%%%%%%%%% bleach correction

%bleachCorrGreen = bleachCorrGreen(1:p);
%bleachCorrRed = bleachCorrRed(1:p);

greenStack = double(greenStack);
redStack = double(redStack);

%rotateDegree = 90;

%greenStack = imrotate(greenStack,rotateDegree);
%redStack = imrotate(redStack,rotateDegree);


%imshow(B(:,:,15),[])

test  = redStack(:,:,1);
figure(1)
imshow(test(:,:,:),[]);
%marker = testNorm;
testNorm = normalizedImage(test);
marker = false(size(testNorm));


for lauf = 1:2

fh = figure(1);
 title('Mark the Background')
frapBkg = roipoly(testNorm(:,:,1));
close(fh);
                

                marker(frapBkg) = 255;
                
end

testNorm = normalizedImage(redStack(:,:,15));


h = fspecial('gaussian', 5, 5) 
redStackGauss = imfilter(redStack(:,:,:),h);
%imshow(redStackGauss(:,:,1),[])
[m n p]= size(redStack);

for lauf = 1:p
             testNorm = normalizedImage(redStackGauss(:,:,lauf));

                mp = imimposemin(testNorm,marker);
                L2 = watershed(mp);

                %mpReg = imregionalmin(mp);
                t3 = false(size(testNorm));
                t3(L2 == 0 ) = 255;
                t3Store(:,:,lauf) = t3;
               % imshow(mpReg(:,:,1),[]);
                
end
                
              %tiffwrite_mat(t3Store, [folderName,'_SegmentGauss.tif']);
    

