function [rgbI] = pol2color(aniso,orient,avg,cmap,varargin)
% Map polarization information to RGB color image.
% rgbI = pol2color(aniso,orient,avg,cmap), where cmap can be 'hsv' or
% 'sbm'.
% rgbI = pol2color(...,'legend',true) inserts a legend in the top-left
% corner.
% rgbI = pol2color(...,'border',pixels) avoids the border pixels in
% normalizing the anisotropy/average intensity.
% 
% Note: Ceiling needs to be applied before calling this function.


% Useful guide on colorspaces:
% http://www.mathworks.com/matlabcentral/fileexchange/28790-colorspace-transformations/content/colorspace/colorspace.html

% The convention I find the most intuitive is as follows:
% orientation is hue
% saturation is anisotropy.
% brightness is intensity
% Such a mapping is provided by HSV and modified LCH space. 

arg.legend=false;
arg.location='top-left';
arg.border=15; % If avgCeiling and anisoCeiling are set to NaN, this value is used to compute the ceilings.
arg.anisoCeiling=NaN; %If anisoCeiling is NaN, set it to max of the data within border.
arg.avgCeiling=NaN; %If avgCeiling is NaN, set it to max of the data within border.
arg=parsepropval(arg,varargin{:});

border=round(arg.border);


if isnan(arg.avgCeiling) 
    avgcrop=avg(border:end-border,border:end-border,:);
    avgCeiling=max(avgcrop(:));
else
    avgCeiling=arg.avgCeiling;
end

if isnan(arg.anisoCeiling)
    anisocrop=aniso(border:end-border,border:end-border,:);
    anisoCeiling=max(anisocrop(:));   
else
   anisoCeiling=arg.anisoCeiling;
end


% Apply ceilings and normalize to 1, for coversion to colorspaces.
aniso=aniso/anisoCeiling;
avg=avg/avgCeiling;
aniso(aniso>1)=1;
avg(avg>1)=1;
orient=orient/pi;
% Incoming matrices may be 3D stacks or movies. Reshape them as 2D
% images for use with colorspace function.
[Y,X,Z]=size(aniso);
% Reshape and scale between 0 and 1 and convert them into image.
% Before returning, reshape rgbI as X*Y*3*Z matrix.
aniso=reshape(aniso,Y,[],1);
orient=reshape(orient,Y,[],1);
avg=reshape(avg,Y,[],1);


switch(cmap)
    case 'hsv' %hue-saturation-value
        hsvI=cat(3,orient,aniso,avg); %hue-saturation-intensity : natural to interprete for anisotropy data and used previously.
        rgbI=hsv2rgb(hsvI);
       % hsvI2=colorspace('rgb->hsv',rgbI);
       % rgbI=colorspace('hsv->',hsvI);
    case 'hs'  %hue-saturation: orientation and anisotropy (set intensity to 1).
        hsI=cat(3,orient,aniso,ones(size(avg)));
        rgbI=hsv2rgb(hsI);
        %rgbI=colorspace('hsv->rgb',hsI);
    case 'sv'  %saturation-intensity: anisotropy and average (color is set to green).
        svI=cat(3,0.33*ones(size(orient)),aniso,avg);
        rgbI=hsv2rgb(svI);
       % rgbI=colorspace('hsv->rgb',siI);
    case 'hv' % orientaiton and average.
        hiI=cat(3,orient,ones(size(aniso)),avg);
        rgbI=hsv2rgb(hiI);
       % rgbI=colorspace('hsv->rgb',hiI);
    case 'hsi' 
        hsiI=cat(3,orient,aniso,avg);
        rgbI=colorspace('hsi->rgb',hsiI);
    case 'lch' % Most natural to interprete.
        % We want highest luminescence with color to imply saturated and
        % oriented pixel. In Lch colorspace, the highest luminescence looks
        % white, irrespective of the color. Therefore, we clip the mapping
        % of luminescence between 0 to 0.5. For making the gray look white,
        % we multiple the resulting rbb image by 2. In effect, we use a
        % colorspace that is a hemisphere of Lch colorspace scaled along
        % the luminescence dimension by factor 2.
        lchI=cat(3,avg,aniso,orient); 
        rgbI=2*colorspace('lch->rgb',lchI);
    case 'hsvG' %Vertical is green, horizontal is magenta. Green is 0.33 in default hsv map.
        orient=mod(orient-0.5+0.33,1); % Turn the wheel by 90 degrees and set the zero to green.
        hsvI=cat(3,orient,aniso,avg); %hue-saturation-intensity : natural to interprete for anisotropy data and used previously.
        rgbI=hsv2rgb(hsvI);
    case 'sbm'
        polI=cat(3,orient,aniso,avg);
        rgbI=pol2opposites(polI);
    otherwise
        error(['pol2color conversion scheme: ' cmap ' is not available.']);
end
        
rgbI=reshape(rgbI,[Y X Z 3]);% Bring the stack back to its original shape.
rgbI=permute(rgbI,[1 2 4 3]);% Permute to make it Y*X*RGB*Z stack, which can be exported using saveastiff.

% rgbI=permute(rgbI,[1 2 4 3]);  % Reshape to make a Y*X*Z*RGB stack. Which vol3d can use for display.

% Color legend for HSV maps.
    if(arg.legend && ~isempty(strfind(cmap,'h'))) % If hue is part of the colormap.
        [x,y]=meshgrid(-1:0.05:1,-1:0.05:1); 
        [theta,rho]=cart2pol(x,-y); % Need to flip Y axis because Y axis runs top to bottom in images.
        theta=mod(theta,pi);
        theta=theta/pi;  %Normalize to 1.

        
        bright=0.9*double(rho<1)+0.1;
        rho(rho>=1)=0;
        legendPol=cat(3,theta,rho,bright);
        legend=hsv2rgb(legendPol);
        
        % Generate the legend using the same scheme as data.
        % Resize legend to 10% of the image.
        legend=imresize(legend,[floor(0.1*X) floor(0.1*Y)]);
        [Ylegend,Xlegend,~]=size(legend);

        % Put the legend in top-left of all slices. Legend is itself an RGB
        % image.
        switch(arg.location)
            case 'top-left'
                rgbI(1:Ylegend,1:Xlegend,:,:)=repmat(legend,[1 1 1 Z]);
            case 'bottom-right'
                rgbI(end-Ylegend+1:end,end-Xlegend+1:end,:,:)=repmat(legend,[1 1 1 Z]);
            case 'top-right'
                rgbI(1:Ylegend,end-Xlegend+1:end,:,:)=repmat(legend,[1 1 1 Z]);
            case 'bottom-left'
                rgbI(end-Ylegend+1:end,1:Xlegend,:,:)=repmat(legend,[1 1 1 Z]);
        end
    end
% Color legend for my custom map that uses opposite colors.    
    if(arg.legend && ~isempty(strfind(cmap,'sbm'))) % If hue is part of the colormap.
        [x,y]=meshgrid(-1:0.05:1,-1:0.05:1); 
        [theta,rho]=cart2pol(x,-y); % Need to flip Y axis because Y axis runs top to bottom in images.
        theta=mod(theta,pi);
        theta=theta/pi;  %Normalize to 1.

        
        bright=0.8*double(rho<1)+0.1;
        rho(rho>=1)=0;
        legendPol=cat(3,theta,rho,bright);
        legend=pol2opposites(legendPol);
        
        % Resize legend to 10% of the image.
        legendsize=floor(0.1*min(X,Y));
        legend=imresize(legend,[legendsize legendsize]);
        [Ylegend,Xlegend,~]=size(legend);

        % Put the legend in top-left of all slices. Legend is itself an RGB
        % image.
        rgbI(1:Ylegend,1:Xlegend,:,:)=repmat(legend,[1 1 1 Z]);
    end
    
    % Since HSV and RGB colormaps are not exactly equivalent, RGB values
    % may be beyond 0 and 1, which can give rise to display errors.
    rgbI(rgbI<0)=0;
    rgbI(rgbI>1)=1;
end