function  plotfxyflip(f)
%The composition of the functions plotfyflip and plotfxflip

xpoints=length(f(:, 1:1))/3;
ypoints=length(f(1:1, :)');
X1=zeros(xpoints,2*ypoints-1); Y1=zeros(xpoints,2*ypoints-1); Z1=zeros(xpoints,2*ypoints-1); 
    X1(:, ypoints:2*ypoints-1)= f(1:xpoints,:);    
    Y1(:, ypoints:2*ypoints-1)= f(xpoints+1:2*xpoints,:); 
    Z1(:, ypoints:2*ypoints-1)= f(2*xpoints+1:3*xpoints,:); 
for j=1:ypoints    
    X1(:, j:j)= f(1:xpoints, ypoints-j+1:ypoints-j+1);
    Y1(:, j:j)= -f(xpoints+1:2*xpoints, ypoints-j+1:ypoints+1-j);
    Z1(:, j:j)= f(2*xpoints+1:3*xpoints, ypoints-j+1:ypoints-j+1);
end

X=zeros(2*xpoints-1,2*ypoints-1); Y=zeros(2*xpoints-1,2*ypoints-1); Z=zeros(2*xpoints-1,2*ypoints-1);

    X(xpoints:2*xpoints-1, :)= X1(1:xpoints,:);  
    Y(xpoints:2*xpoints-1, :)= Y1(1:xpoints,:);
    Z(xpoints:2*xpoints-1, :)= Z1(1:xpoints,:);
for j=1:xpoints    
    X(j:j, :)= -X1(xpoints-j+1:xpoints-j+1,:);
    Y(j:j, :)= Y1(xpoints-j+1:xpoints-j+1,:);
    Z(j:j, :)= Z1(xpoints-j+1:xpoints-j+1,:);
end

surf(X,Y,Z,'Clipping','off','EdgeAlpha',0.5,'FaceAlpha',1,'FaceLighting','gouraud','EdgeLighting','gouraud','EdgeColor','white','FaceColor','b')
daspect([1 1 1])
grid off
end