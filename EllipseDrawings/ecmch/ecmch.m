function f = ecmch(pot, IC, Ix, Iy, looporder,H)
%Update: 12 Nov 2014
%Computes a CMC surface via the DPW (generalized Weierstrass) method
%from the potential pot.   The computation is done over a rectangular grid,
%determined by Ix and Iy.  All loops will be expanded up to order
%'minlooporder' (or the maximum order of the other input data in the loop parameter,
%whichever is the greater) on the inner bands, and 'maxlooporder' on the outer bands. 
%For small regions try minlooporder=maxlooper=4.  For larger regions, something like
%minlooporder=4, maxlooporder=8.  Usually the higher loop order is needed only on the outer
%bands, so time is saved by doing this.  Increasing the loop orders will result in greater
%accuracy, but longer computation time.  Also plots the surface using surf.

%The loop group construction, with the same setup used here, is
%described in 	 arXiv:0908.3274  [math.DG].  The loops described there
%have to be untwisted before entering them into this program.  This is
%achieved by: twisted  
%        [a(lambda),  b(lambda);  c(lambda),  d(lambda)]
%           ->
%        [a(sqrt(lambda)),  B_-(sqrt(lambda));  C_+(sqrt(lambda)), d(sqrt(lambda)]   
%  untwisted 
% where B_-(lambda):= lambda^{-1} b(lambda) and
%C_+(lambda):= lambda*c(lambda).


%Author: David Brander;  email: D.Brander@mat.dtu.dk.


% Input: All loops are untwisted. 
% pot = @(z,h) [A_-n(z,h),... A_n(z,h)]
% where each A_i(z,h) is a 2 x 2 matrix valued expression involving z and
%h, and holomorphic in z.
%This represents the potential (A_-n lambda^{-1} + A_1 lambda + ... + A_n
%lambda^n) dz. 

% All loops are converted into large square matrices, such that every pair of rows is 
%a copy of the previous pair, shifted to the right by 2 columns. Then matrix
% multiplication corresponds exactly to loop multiplication, and the LU
% decomposition corresponds to the birkhoff decomposition of loops.

%IC  represenst the initial condition, and is a Laurent
%polynomial (untwisted) loop: IC = [F_-j , ...., F_j]
%    where F_i  are constant 2x2  matrices.  Usually IC = eye(2).

%Ix should be of the form [x_0, stepsize, n1, n2],  where n1 is the number
%of steps to the left of x_0 and n2, the number to the right.  Similar for
%Iy.

%"looporder" can be a positive integer "looporder=maxlooporder"
%or a pair "looporder=[maxlooporder, minlooporder]" 
%(Default is minlooporder = floor(maxlooporder/4)).  Loops
%are approximated as polynomials of order minlooporder at the initial 
%point up to maxlooporder at the points furthest from the initial point.
%Default for minlooporder is the 1 or the maximum order of Fp0 and Fm0.

% H is the mean curvature.  The Sym-formula is divided by H.
%the potential pot should be a function handle @(z,h), where h is the mean
%curvature.

%OUTPUT:
%f is of the form [X; Y; Z], where each of X,Y and Z are m x n matrices,
%m is the number of x gridpoints 1+Ix(3)+Ix(4), and n the number of y
%grid points 1+Iy(3)+Iy(4).
%X is the matrix of x-values at the grid
%points, Y the matrix of y values, and Z the matrix of z values. 
%This can be plotted with the command 
%   surf(f(1:m, :), f(m+1: 2*m, :),  f(2*m+1: 3*m, :)),
% (that is, surf(X,Y,Z) ) which is also done automatically.

%Error estimates are computed by checking if the solutions matrices are in
%su(2).  The formula used is norm(f'+f)/norm(f),  where f'+f=0 is the su(2)
%condition. The error terms should be less than about 10^(-2), otherwise, a 
%higher looporder should be used.

lambda=1;
A = @(z)pot(z,H);
%checkinput parses the input for consistency, converts the loops into
%square matrix form and sets the minimum and maximum looporders.
[xpoints,ypoints,bigA,bigIC,minlooporder,maxlooporder,m]=checkinput(A, IC, Ix , Iy, looporder);
looporders = spreadlooporders(Ix,Iy,minlooporder,maxlooporder); %the orders for the LU decomposition and loop evaluation
bsize=m*(2*maxlooporder+1);  %dimensions of the square matrix form of loops
%first integrate along the line y=0:
F1=realSintegrateA(bigA, bigIC, Ix,m);
%then integrate along each line x=x_i:
xline = makexline(Ix);
Y=zeros(2*xpoints,2*ypoints);
ErrorMatrix=zeros(xpoints,ypoints);
for row=1:xpoints
    Phi=imagSintegrateA(bigA,bigIC, Iy, xline(row),m);
    for j=1:ypoints
        n=looporders(row,j);
       k=min([maxlooporder,n]);
       Phimat=loop2SMat(Phi(:, (j-1)*bsize+1:j*bsize),maxlooporder,k);
       thisPhiloop=F1(:, (row-1)*bsize+1:row*bsize)*Phimat;         
       thisPhi=loop2SMat(thisPhiloop,n,k);
       [V,~]=lu2(thisPhi'*thisPhi,m,0.00000001,n);
       U=loop2SMat(V,maxlooporder,n);
       F=thisPhiloop/U;
       [fval, UError] = sym2(F,H,lambda);    
       Y(m*row-1:m*row,m*(j-1)+1:m*j)=fval;
       ErrorMatrix(row,j)=UError;
    end
    if mod(row,10)==1
        disp(['Row ', num2str(row),  ' Max Error ', num2str(max(ErrorMatrix(row:row,:)),2),'. Errors: ', num2str(ErrorMatrix(row:row,:),2)]); 
    end
  clear Phi; 
end
 disp(' ');
disp(['Max error:', num2str(max(max(ErrorMatrix)),2), '. Mean error: ', num2str(mean(mean(ErrorMatrix)),2), '. Corner Errors: ']);
  disp(num2str([ErrorMatrix(1,1), ErrorMatrix(1,ypoints);ErrorMatrix(xpoints,1), ErrorMatrix(xpoints,ypoints)],2));
 f = surfplotdata(Y);
end


     %%%%%%%%%%%%%%%%%%%%%%%%%%%% sub functions %%%%%%%%%%%%%%%%
     

           
            function [V,p]=lu2(X,m,thresh,n)
   %computes anlu decomposition of X,
   %where X is the small square matrix form of a loop.
   %P is 0 if X is in the big cell
[~,Umat,P]=lu(X(1:m*(n+1),1:m*(2*n+1)),thresh);
p=trace(P)-length(P);
  UmDims=size(Umat);
  V=Umat(UmDims(1)-m+1:UmDims(1),1:m*(2*n+1));
           end
 
      
    function xline = makexline(Ix) 
     %makes the interval of integration along the x axis,
xline=zeros(1,1+Ix(3)+Ix(4));
xline(Ix(3)+1) = Ix(1);
for count=1:Ix(3)
    xline(count)=Ix(1)-Ix(2)*(Ix(3)-count+1);
end
for count=1:Ix(4)
    xline(Ix(3)+1+count)=Ix(1)+Ix(2)*count;
end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [xpoints,ypoints,bigA,bigIC,minlooporder,maxlooporder,m]=checkinput(A, IC, Ix , Iy, looporder);
%function [xpoints,ypoints,bigAm, bigAp, bigFm0, bigFp0,minlooporder,maxlooporder]=checkinput(Am, Ap, Fm0, Fp0, Iy , Ix, looporder)

%Checks the input data for consistency with expected input.
%Converts the intial data, in terms of simple untwisted loops, into 
%untwisted loops represented by big square matrices, expanded up to order
%'looporder' or the maximum order of the other data in the loop parameter,
%whichever is the greater.

% Output:  Each of bigAp and bigAm is a handle for a function of the form  = @(t) M(t)
% where M is a  2(4n+1) x 2(4n+1) square matrix valued function of t
% which looks something like 
%      [x_0 ... x_n 0 0 ... 0;
%       x_-1 x_0 ....x_n 0 ...;
%       x_-2 x_-1 x_0 .......; ..... etc,  each x_i a 2x2 matrix

% bigFp0 and bigFm0 are square matrices of initial values, of the same size 
%and similar structure to A.
assert((Ix(3)>=0)&(Ix(4)>=0)&(Ix(3)+Ix(4)>0),'Ix(3) and Ix(4) must be non-negative, and at least one positive');
assert((Iy(3)>=0)&(Iy(4)>=0)&(Iy(3)+Iy(4)>0),'Iy(3) and Iy(4) must be non-negative, and at least one positive');
xpoints=1+Ix(3)+Ix(4); ypoints=1+Iy(3)+Iy(4);

ICDims=size(IC); 
ADims=size(A(Ix(1)+1i*Iy(1)));
m=ICDims(1);
assert((ICDims(1)==2)&(ADims(1)==2),'All input loops must be matrices with two rows');

ICorder = (ICDims(2)-m)/(2*m);  
AOrder = (ADims(2)-m)/(2*m);  
assert(norm(floor([ICorder, AOrder])-[ICorder, AOrder])==0,...
    'Number of columns for input matrices and function handles are not correct for loops');
loDims=size(looporder);
assert((loDims(2)<=2)&(norm(floor(looporder)-looporder)==0)&(norm(looporder-abs(looporder))==0),...
    'looporder must be a non-negative integer or a 1x2 non-negative integer matrix ');
if loDims(2)==2
    maxlooporder=max([looporder(1),ICorder,AOrder]);
    minlooporder=looporder(2);
else
  maxlooporder=max([looporder,ICorder,AOrder]); 
  minlooporder=max([floor(maxlooporder/4),ICorder]);  
end
%convert the loops into square matrix form:
bigA=@(t)loop2SMat(A(t),maxlooporder,AOrder);
bigIC= loop2SMat(IC,maxlooporder,ICorder);
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%

    function orders = spreadlooporders(Ix,Iy,minlooporder,maxlooporder)
 %creates a matrix of looporders to be used in the LU decomposition
 %to speed things up.
 
orders=zeros(Ix(3)+Ix(4)+1,Iy(3)+Iy(4)+1);
center=[Ix(3)+1,Iy(3)+1];            %the coordinates of the center of the matrix
cornerdistances=[norm(center-[1,1]), norm(center-[1,Iy(3)+Iy(4)+1]), norm(center-[Ix(3)+Ix(4)+1,1]),  norm(center-[Ix(3)+Ix(4)+1,Iy(3)+Iy(4)+1])];
for m=1:(Ix(3)+Ix(4)+1)
   for k=1:Iy(3)+Iy(4)+1
       orders(m,k)=min(ceil(minlooporder+(maxlooporder-minlooporder)*norm(center-[m,k])/max(cornerdistances)),maxlooporder);
   end
end
     end 
%%%%%%%%%%%%%%%%%%%%%%%%%%
   
     
    function F = realSintegrateA(A, F0, Ix,m)
%solves the equation dF/dt = F A, with initial condition F(t0) = F0, with t
%in the specified interval
%The values of F are returned in loops in row form

tempIx3=Ix(3); tempIx4=Ix(4);
if Ix(3)<2, tempIx3=2; end
if Ix(4)<2, tempIx4=2; end
stepy= Ix(2);  intervalyl=zeros(1, tempIx3+1); intervalyr=zeros(1,tempIx4+1); 
for count = 1 : tempIx4+1
    intervalyr(count)=  Ix(1) + (count-1)*stepy;
end
for count=1: tempIx3+1
    intervalyl(count)=  Ix(1) - (count-1)*stepy;
end
bsize=length(F0);
n =(bsize/m-1)/2; % the order of X in lambda
initF=reshape(F0(m*n+1:m*n+m,:),m*bsize,1);
testfn = @(t,y) reshape(sparse(reshape(y,m,bsize))*A(t), m*bsize,1);
[~,soll] = ode45(testfn ,intervalyl,initF); 
[~,solr] = ode45(testfn ,intervalyr,initF);  
F=sparse([reshape(flipud(soll(2:Ix(3)+1, :)).', m,bsize*(Ix(3))),   reshape(solr(1:Ix(4)+1, :).', m,bsize*(Ix(4)+1))]);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function F = imagSintegrateA(A, F0, Iy, xvalue,m)
%The integration along the y axis.
%solves the equation dF/dt = F A(xvalue + i*t), with initial condition F(t0) = F0, with t
%in the specified interval
%The values of F are returned in loops in row form

tempIy3=Iy(3); tempIy4=Iy(4);
if Iy(3)<2, tempIy3=2; end
if Iy(4)<2, tempIy4=2; end
stepy= Iy(2);   intervalyr=zeros(1,tempIy4+1); intervalyl=zeros(1,tempIy3+1);
for j = 1 : tempIy4+1
    intervalyr(j)=  Iy(1) + (j-1)*stepy;
end
for j = 1 : tempIy3+1
    intervalyl(j)=  Iy(1) - (j-1)*stepy;
end
bsize=length(F0);
n =(bsize/m-1)/2; % the order of X in lambda
initF=reshape(F0(m*n+1:m*n+m,:),m*bsize,1);
testfn = @(t,y) reshape(sparse(reshape(y,m,bsize))*(1i*A(xvalue+1i*t)), m*bsize,1);
if Iy(3)>0
  [~,soll] = ode45(testfn ,intervalyl,initF);
  F1=sparse(reshape(flipud(soll(2:Iy(3)+1, :)).', m,bsize*(Iy(3))));
else
    F1=[]; 
end
if Iy(4)>0
  [~,solr] = ode45(testfn ,intervalyr,initF);
  F2=sparse(reshape(solr(2:Iy(4)+1, :).', m,bsize*(Iy(4))));
else
    F2=[]; 
end
F=sparse([F1, F0(m*n+1:m*n+m,:) F2]);
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,UError] = sym2(X,H,lambda)
%computes the the Sym formula 
%S(X,r) = -2*i*r^2 Ad_D (d/d lambda X)(r^2) X(r^2)^{-1}
%at lambda =r=1 , where X should be the untwisted frame.

%D is the matrix [sqrt(lambda), 0;  0, sqrt(lambda)^(-1)]
%Input: X is a block diagonal matrix X=[B_1,.......,B_points], where each block
%B_i is itself a 2x2  block diagonal B_i = [s_i, s_i, .. s_i], where s_i is
%the value of f at point i.
n =(length(X)-2)/4; % the order of X in lambda
Xdiff=zeros(2,2);
  for j=1:n
      Xdiff=Xdiff+j*lambda^(j-1)*X(:,2*(n+j)+1:2*(n+j+1));
  end
  for j=1:n
      Xdiff=Xdiff-j*lambda^(-j-1)*X(:,2*(n-j)+1:2*(n-j+1));
  end
Xval=rloopeval(X,lambda);
f=-(2*1i/H)*Xdiff/Xval;
UError=0;
if norm(f)>0
  UError=norm(f'+f)/norm(f); %check to see if the matrix is in U(2)
end
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = surfplotdata(f)
%converts the solution f into a standard form (X;Y;Z)
%for plotting with the surf command
fheight=length(f(:,1:1));
xpoints=fheight/2;

fwidth=length(f(1:1,:));
ypoints=fwidth/2;

Y=zeros(3*xpoints, ypoints);
for k=1:xpoints
    for j=1:ypoints  
        Y(k,j) = -imag(f(2*k-1, 2*j));      
        Y(xpoints + k, j) = real(f(2*k-1, 2*j));
        Y(2*xpoints+k, j) =imag(f(2*k-1, 2*j-1));  
    end
end
        
surf(Y(1:xpoints, :), Y(xpoints+1: 2*xpoints, :),  Y(2*xpoints+1: 3*xpoints, :),'Clipping','off','EdgeAlpha',0.4,'FaceAlpha',1,'FaceLighting','gouraud','EdgeLighting','gouraud','DiffuseStrength',0.8,'EdgeColor','white','FaceColor',[0.05,0.55,0.85])
 daspect([1 1 1])
 grid off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = rloopeval(X,r)
% Evaluates a loop at lambda = r.
% Input: X = [x_-n .... x_0 x_1 .... x_n], where each x_i is 2x2
%  r is a number.
%Output: a sparse 2x2 matrix, the sum of x_i r^i
n = (length(X)-2)/4; % the order of X in lambda

Y= X(1:2, 2*n+1:2*n+2);
for count = 1 : n 
    Y = Y +  r^count*X(1:2, 2*n+1+2*count : 2*n+2+ 2*count) +  r^(-count)*X(1:2, 2*n+1-2*count : 2*n + 2 - 2*count);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function Xmat = loop2SMat(X,n,k)
%Converts a loop into a big square matrix 
%Input:  a matrix X=[x_-nold .... x_0 x_1 .... x_nold] representing a Laurent
%polynomial loop of order nold
%The input loop does not need to be the same order as the output loop

%x_-nold lambda^-nold + ... + x_0 + x_nold lambda^nold, where each x_i is an mxm matrix 
%Output:  Xmat is a sparse m(4n+1) x m(4n+1) square matrix, which looks something
%like  [x_0 ... x_n 0 0 ... 0;
%       x_-1 x_0 ....x_n 0 ...;
%       x_-2 x_-1 x_0 .......; ..... etc
% Multiplying  Xmat with Ymat corresponds to multiplying the loops X and Y.
Xdims=size(X);
m=Xdims(1);
nold =(Xdims(2)/m-1)/2; % the order of the input loop X in lambda
k=min(nold,k);
Xmat=zeros(m*(2*n+1));
Xmat(m*n+1:m*(n+1),m*(n-k)+1:m*(n+1+k))=full(X(1:m, m*(nold-k)+1: m*(nold+k+1)));
for j=0:n-1
    Xmat(m*j+1:m*j+m,1:m*(n+1+j))=Xmat(m*n+1:m*(n+1),m*(n-j)+1:m*(2*n+1));
    Xmat(m*(2*n-j)+1:m*(2*n-j+1),m*(n-j)+1:m*(2*n+1))=Xmat(m*n+1:m*(n+1),1:m*(n+1+j));
end
Xmat=sparse(Xmat);
end




