function updatePolHist(p,ha,hfig, aniso, orient, avg, I0, I135, I90, I45)
persistent  hLineROI hLineOrient;
% p is ROI, ha is handles to the axes.
    histMask=poly2mask(p(:,1),p(:,2),size(I0,2),size(I0,1));
   rpMask=regionprops(histMask,'Orientation','MajorAxisLength','Centroid','MinorAxisLength');
   Lm=rpMask.MajorAxisLength; % Length of line indicating orientation of the mask.
   maskCen=rpMask.Centroid;
   maskOrient=mod(rpMask.Orientation*(pi/180),pi);
    % Determine the display limits after ignoring the
    % border.
    I0crop=I0(histMask);
    I135crop=I135(histMask);
    I90crop=I90(histMask);
    I45crop=I45(histMask);
    OrientationCrop=orient(histMask);
    AnisotropyCrop=aniso(histMask);
    AverageCrop=avg(histMask);
    maxAniso=max(AnisotropyCrop(:));
    Ivec=[I0crop(:); I135crop(:); I90crop(:); I45crop(:)];
    Ilims=[min(Ivec) max(Ivec)];
    
   [meanOrient,meanAniso]=ComputeFluorAnisotropy(mean(I0crop),mean(I45crop),mean(I90crop),mean(I135crop),'anisotropy');
   Lf=Lm*meanAniso; % Scale down the lenght of line indicating orientation of fluorophore depending on the anisotropy.
   
    set([ha(1:4) ha(7)],'Clim',Ilims);
    set(ha(6),'Clim',[0 maxAniso]);
    axes(ha(6)); title(['Anisotropy: max ' num2str(maxAniso,2) ]);
 
    axes(ha(7));

    
    % Note that order of Y-coordinates is flipped because the axes follow
    % the image convention (positive from top to bottom).
    % Draw a line indicating the ROI.
    if(ishandle(hLineROI))
        delete(hLineROI);
    end
    
    hLineROI=line([maskCen(1)-0.5*Lm*cos(maskOrient) maskCen(1)+0.5*Lm*cos(maskOrient)],...
        [maskCen(2)+0.5*Lm*sin(maskOrient) maskCen(2)-0.5*Lm*sin(maskOrient)],'LineWidth',2);
    
    % Draw a line indicating the polarization orientation.
    
    if(ishandle(hLineOrient))
        delete(hLineOrient);
    end
    hLineOrient=line([maskCen(1)-0.5*Lf*cos(meanOrient) maskCen(1)+0.5*Lf*cos(meanOrient)],...
        [maskCen(2)+0.5*Lf*sin(meanOrient) maskCen(2)-0.5*Lf*sin(meanOrient)],'LineWidth',2,'Color','Red');
    
   
    subplot(2,4,8,'Parent',hfig); 
    polarPlotAnisoStat(AnisotropyCrop,OrientationCrop,AverageCrop,'Nbins',30,'PlotType','Polar','ReferenceOrient',maskOrient);
    TotalI=mean(AverageCrop(:));
    SNR=sqrt(TotalI);
    title(['Intensity vs. Orientation (mean aniso:' num2str(meanAniso,2) ', SNR:' num2str(SNR,3) ')'] ); 
%     [counts,levels]=hist(mod(OrientwrtMask,180),0:0.5:180);
%     stem(levels,counts);
%    
%     xlim([0 180]); 
%     

end    
