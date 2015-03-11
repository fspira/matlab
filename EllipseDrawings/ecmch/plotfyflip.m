function  plotfyflip(f)
%for surfaces which are symmetric with respect to the xz plane.
%If f is a half of this surface, lying to one side of this plane,
%the output is the whole surface. 
%Input:  f is assumed to be a matrix of
%the form (X;Y;Z), as in the input for the function surf.

xpoints=length(f(:, 1:1))/3;
ypoints=length(f(1:1, :)');
X=zeros(xpoints,2*ypoints-1); Y=zeros(xpoints,2*ypoints-1); Z=zeros(xpoints,2*ypoints-1);


 
    
    
    X(:, ypoints:2*ypoints-1)= f(1:xpoints,:);    
    Y(:, ypoints:2*ypoints-1)= f(xpoints+1:2*xpoints,:); 
    Z(:, ypoints:2*ypoints-1)= f(2*xpoints+1:3*xpoints,:); 
for j=1:ypoints    
    X(:, j:j)= f(1:xpoints, ypoints-j+1:ypoints-j+1);
    Y(:, j:j)= -f(xpoints+1:2*xpoints, ypoints-j+1:ypoints+1-j);
    Z(:, j:j)= f(2*xpoints+1:3*xpoints, ypoints-j+1:ypoints-j+1);
end
surf(X,Y,Z,'Clipping','off','EdgeAlpha',0.5,'FaceAlpha',1,'FaceLighting','gouraud','EdgeLighting','gouraud','EdgeColor','white','FaceColor','b')
daspect([1 1 1])
grid off
end