function f= ecmchpolar(pot, IC, Ir, It, looporder,H)
%Updated: 13 Nov 2014
% Similar to ecmch, only integrates over a polar region, starting from
%the center. The potntial,
% pot, is given as a function handle of (r,t), instead of z, where z=r
% exp(it).  

%Ir should be of the form [0, stepsize, n],  where n is the integration is
%over the interval [0, 0+ stepsize*(n+1)].
% It is of the form [t_0, stepsize, m], where t_0 is the initial angle, 
% and the integration is over [t_0,  stepsize*(m+1)]

%H is the mean curvature.  The Sym-formula is divided by H.
%the potential pot should be a function handle @(r,t,h), where h is the mean
%curvature.
lambda=1;  %can replace with any unit length complex number
A = @(r,t)pot(r,t,H);
%checkinput parses the input for consistency, converts the loops into
%square matrix form and sets the minimum and maximum looporders.
[rpoints,tpoints,bigA,bigIC,minlooporder,maxlooporder,m]=checkinput(A, IC, Ir , It, looporder);
looporders=setlooporders(Ir,minlooporder,maxlooporder);
bsize=m*(2*maxlooporder+1);  %dimensions of the square matrix form of loops
% integrate along the line t=t0+j*stepsize:
Y=zeros(2*rpoints,2*tpoints);
ErrorMatrix=zeros(tpoints,rpoints);
for row=1:tpoints
   Phi=realSintegrateA(bigA, bigIC, Ir, It(1)+(row-1)*It(2),m);
    for j=1:rpoints
       n=looporders(j);
       k=min([maxlooporder,n]);
       Phimat=loop2SMat(Phi(:, (j-1)*bsize+1:j*bsize),n,k);
       thisPhiloop=bigIC(m*maxlooporder+1:m*maxlooporder+m,m*(maxlooporder-n)+1:m*(maxlooporder+1+n))*Phimat;
       thisPhi=loop2SMat(thisPhiloop,n,k);
       [V,~]=lu2(thisPhi'*thisPhi,m,0.00000001,n);
       U=loop2SMat(V,n,n);
       F=thisPhiloop/U;
       [fval, UError] = sym2(F,H,lambda);    
       Y(m*(j-1)+1:m*j,m*row-1:m*row)=fval;
       ErrorMatrix(row,j)=UError;
    end
    if mod(row,10)==1
        disp(['Row ', num2str(row),  ' Max Error ', num2str(max(ErrorMatrix(row:row,:)),2),'. Errors: ', num2str(ErrorMatrix(row:row,:),2)]); 
    end
  clear Phi; 
end
 disp(' ');
disp(['Max error:', num2str(max(max(ErrorMatrix)),2), '. Mean error: ', num2str(mean(mean(ErrorMatrix)),2), '. Corner Errors: ']);
  disp(num2str([ErrorMatrix(1,1), ErrorMatrix(1,rpoints);ErrorMatrix(tpoints,1), ErrorMatrix(tpoints,rpoints)],2));
 f = surfplotdata(Y);
end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     %%%%%%%%%%%%%%%%%%%%%%%%%%%% sub functions %%%%%%%%%%%%%%%%

            function [V,p]=lu2(X,m,thresh,n)
   %computes an lu decomposition of X,
   %where X is the small square matrix form of a loop.
   %P is 0 if X is in the big cell
[~,Umat,P]=lu(X(1:m*(n+1),1:m*(2*n+1)),thresh);
p=trace(P)-length(P);
  UmDims=size(Umat);
  V=Umat(UmDims(1)-m+1:UmDims(1),1:m*(2*n+1));
           end
 



    function F = realSintegrateA(A, F0, Ix,t0,m)
    %solves the equation dF/dt = F A, with initial condition F(t0) = F0, with t
%in the specified interval
%The values of F are returned in loops in row form

tempIx3=max(Ix(3),2); 
interval=0:Ix(2):Ix(2)*(tempIx3);

bsize=length(F0);
n =(bsize/m-1)/2; % the order of X in lambda
initF=reshape(F0(m*n+1:m*n+m,:),m*bsize,1);
testfn = @(r,y) reshape(sparse(reshape(y,m,bsize))*exp(1i*t0)*A(r,t0), m*bsize,1);
[~,sol] = ode45(testfn ,interval,initF); 
F=sparse(reshape(sol(1:Ix(3)+1, :).', m,bsize*(Ix(3)+1)));  

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
%n =(length(Xmat)-2)/8; % the order of X in lambda
%X= Xmat(length(Xmat)/2: length(Xmat)/2 +1,  2*n+1 : 6*n+2); %convert matrix to loop 
n=(length(X)-2)/4;
Xdiff=zeros(2,2);
  for j=1:n
      Xdiff=Xdiff+j*lambda^(j-1)*X(:,2*(n+j)+1:2*(n+j+1));
  end
  for j=1:n
      Xdiff=Xdiff-j*lambda^(-j-1)*X(:,2*(n-j)+1:2*(n-j+1));
  end
f=-(2*1i/H)*Xdiff/rloopeval(X,lambda);
UError=0;
if norm(f)>0
  UError=norm(f'+f)/norm(f); %check to see if the matrix is in U(2)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = surfplotdata(f)
%
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



function Y = rloopeval(X,r)
% Evaluates a loop at lambda = r.
% Input: X = [x_-n .... x_0 x_1 .... x_n], where each x_i is 2x2
%  r is a number.
%Output: a sparse 2x2 matrix, the sum of x_i r^i
n = (length(X)-2)/4; % the order of X in lambda
Y= full(X(1:2, 2*n+1:2*n+2));
for count = 1 : n 
    Y = Y +  r^count*X(1:2, 2*n+1+2*count : 2*n+2+ 2*count) +  r^(-count)*X(1:2, 2*n+1-2*count : 2*n + 2 - 2*count);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function looporders=setlooporders(Ir,minlooporder,maxlooporder)

looporders=maxlooporder*ones(1,Ir(3)+1);
if maxlooporder>minlooporder
  rstepsize=(maxlooporder-minlooporder)/Ir(3);
  looporders=min(looporders,ceil((minlooporder/rstepsize:1:maxlooporder/rstepsize)*rstepsize));  
end
end





function [rpoints,tpoints,bigA,bigIC,minlooporder,maxlooporder,m]=checkinput(A, IC, Ir , It, looporder)

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
assert((Ir(3)>0),'Ir(3) must be positive');
assert((It(3)>0),'It(3) must be positive');
rpoints=1+Ir(3); tpoints=1+It(3);

ICDims=size(IC); 
ADims=size(A(Ir(1),It(1)));
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
bigA=@(r,t)loop2SMat(A(r,t),maxlooporder,AOrder);
bigIC= loop2SMat(IC,maxlooporder,ICorder);
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%