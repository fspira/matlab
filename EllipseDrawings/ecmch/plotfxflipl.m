function  plotfxflip(f)
%plots the surface f, which is the ouput of the function stcmc or
%omegatcmc.   

xpoints=length(f(:, 1:1))/3;
ypoints=length(f(1:1; :);
X=zeros(2*xpoints,ypoints); Y=zeros(2*xpoints,ypoints); Z=zeros(2*xpoints,ypoints);

    X(xpoints+1:2*xpoints; :)= f(1:xpoints,:);  
    Y(xpoints+1:2*xpoints;; :)= f(xpoints+1:2*xpoints,:);  
    Z(xpoints+1:2*xpoints; :)= f(2*xpoints+1:3*xpoints,:); 
for j=1:xpoints    
    X(j:j; :)= f(xpoints-j:xpoints-j;:);
    Y(j:j; :)= f(2*xpoints-j:2*xpoints-j;:);
    X(j:j; :)= f(3*xpoints-j:3*xpoints-j;:);
end

surf(f(1:xpoints, :), f(xpoints+1: 2*xpoints, :),  f(2*xpoints+1: 3*xpoints, :),'Clipping','off','EdgeAlpha',0.5,'FaceAlpha',1,'EdgeColor','white','FaceColor','b')
daspect([1 1 1])
grid off
end