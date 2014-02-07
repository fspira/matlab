
clear all

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;

javaaddpath '/Applications/Fiji.app/scripts/Miji.m'
javaaddpath '/Applications/Fiji.app/jars/ij-1.47b.jar'
addpath('/Applications/Fiji.app/scripts')

addpath('/Applications/Fiji.app/scripts')
addpath('/Users/spira/Documents/Matlab_scripte/Image_Processing_utils')
addpath('/Users/spira/Documents/Matlab_scripte/tiffIO')
addpath('/Users/spira/Documents/Matlab_scripte/')
addpath('/Users/spira/Desktop/Desktop/LifeactCherry_GlGPIEgfp/131204')
addpath('/Users/spira/Documents/MATLAB_scripte/ImageProcessing/Utilities')
curdir = pwd;


%%%% If single planes are used, split file into different channels and read
%%%% each channel individually

%tifFilenameGreen = 'cell1_halfStack-grey-P.tif';


%tifFilenameRed = 'cell1_halfStack-grey-S.tif';


%imgOrigGreen = imread(char(tifFilename));


%imgOrigRed = imread(char(tifFilenameRed));



    % truncName = findstr(tifFilename,'.lsm');
    % folderName = tifFilename(1:truncName-1);
    
tifFilename = 'cell12_halfStack.lsm.lsm';
imgOrig = tiffread30(char(tifFilename))
img = cat(3,imgOrig.data);

[m n p] = size(img);

for lauf = 1:p
    imgOrigGreen(:,:,lauf) = img{:,1,lauf};
    imgOrigRed(:,:,lauf) = img{:,2,lauf};
end
%greenImg = imgOrigGreen(:,:,2);
%redImg = imgOrigRed(:,:,1);

greenImg = double(imgOrigGreen);
redImg = double(imgOrigRed);

%[redImg greenImg] = bkgCorrectionRedGreen(redImg, greenImg);

%%%%% Setting negative values to small values, to avoid division through zero
greenNormZero = greenImg > 0;
greenNorm = double((greenImg .* greenNormZero));

redNormZero = redImg > 0;
redNorm = double(redImg .* redNormZero);

%for lauf = 1:p
%    greenNorm(:,:,lauf) = greenNorm(:,:,lauf)./mean(mean(greenNorm(:,:,lauf)));
%    redNorm(:,:,lauf) = redNorm(:,:,lauf)./mean(mean(redNorm(:,:,lauf)));
%end


S1 = greenNorm;
S2 = redNorm;


%%%% Formula according to ZEISS manual G-Factor was calculated using the
%%%% exact same setting and Alexa488 Phalloidin in solution. Ch1 = P =
%%%% green Image, Ch2 = S = red Image

D=((S1-1.079.*S2)./(S1+(2*1.079).*S2)).*255;
%D=((a-1.079.*b)./(a+(2*1.079).*b));

%%% since the background gives rise to artificial 255 values, I remove this
%%% specific value.

%tmp = D< 100;
%D = D .* tmp;


%tmp = D > -80;
%D = D .* tmp;




%%%%%%%%% separate the horizontal from the vertical component
for lauf = 1:p
    minimumValue = min(min(D(:,:,lauf)));
    dTmp = D(:,:,lauf) > minimumValue & D(:,:,lauf) <0;
    dNeg = D(:,:,lauf) .* dTmp;
    dNegInv(:,:,lauf) = -dNeg;
end

%%%%%% Remove negative Values from D
tmp = D >0;
D = D .* tmp;


%%%%%% replace inf with zero
D(~isfinite(D)) = 0;
dNegInv(~isfinite(dNegInv)) = 0;


dNorm = normalizedImage3D(D);
dNegNorm = normalizedImage3D(dNegInv);

for lauf = 1:p
    dNegRange(:,:,lauf) = (dNegInv(:,:,lauf)./max(max(dNegInv(:,:,lauf)))).*127;
    dRange(:,:,lauf) = ((D(:,:,lauf)./ max(max(D(:,:,lauf)))).*126)+128;
    tmp = dRange(:,:,lauf) >129;
    dRange(:,:,lauf) = dRange(:,:,lauf) .*tmp;
end

%ratio = dNorm./dNegNorm;
dSum = dNegRange +dRange;
imshow(dSum(:,:,1),[])
dSum = uint8(round(dSum));
colormap(jet)
figure(1)
imshow(dNorm(:,:,1),[])
figure(2)
imshow(dNegNorm(:,:,1),[])

%imwrite(dSum,[tifFilename,'_Aniso_FalseColor.tif'])
tiffwrite_mat(dSum, [tifFilename,'_Aniso_FalseColor.tif'])

%colormap(jet)
close

tiffwrite_mat(dNegInv, [tifFilename,'_P.tif'])
tiffwrite_mat(D, [tifFilename,'_S.tif'])


for lauf = 1:p

 S1_green= ind2rgb(round(dNorm(:,:,lauf)),green);
 S2_red = ind2rgb(round(dNegNorm(:,:,lauf)),red);
             
    imgMerge(:,:,:,lauf) = S1_green + S2_red;
    
end

tiffwrite_RBG(imgMerge,[tifFilename,'Aniso_RGB.tif'])
%imwrite(imgMerge,[tifFilename,'_Aniso_RGB.tif'])
%imgMergeTmp = imgMerge(:,:,1);
%    imgMergeScale = imgMergeTmp(1:2:end, 1:2:end);

%    meshSize = size(D)
% [x,y] = meshgrid(1:2: meshSize(2),1:2:meshSize(1));

% dScaleTMP = D(:,:,1);
% dNegScaleTmp = -dNegInv(:,:,1);
 
 %       dScale = double(dScaleTMP(1:2:end, 1:2:end));
 %       dNegScale = double(dNegScaleTmp(1:2:end, 1:2:end));
        
        %%%%% normalize to mean intensity
        
  %      dScaleNorm = dScale./mean(mean(dScale));
  %      dNegScaleNorm = -dNegScale./mean(mean(-dNegScale));
   
  %      imshow(imgMergeScale,[]);hold on; quiver(x,y, -dScale(:,:,1),  -dNegScale(:,:,1),4)

%imwrite(D,[tifFilenameGreen,'_Aniso'])

%imwrite(D,jet,[tifFilenameGreen,'_Aniso.tif'],'Compression','none',...
%'WriteMode','append')
close all

    
