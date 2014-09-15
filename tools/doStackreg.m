function [Stack_Reg] = doStackreg(Stack_In)

%Stack_In = GreenBeadPic;

[m n p] = size(Stack_In);

Stack_Reg(:,:,1) = Stack_In(:,:,1);

for subrun = 2 : p
       
        subrun;
     
         [output,Greg] = dftregistration(fft2(Stack_In(:,:,1)),fft2(Stack_In(:,:,subrun)),100);
        row_shift = output(3);
        col_shift = output(4);
        
        
        [nr,nc]=size(Stack_In(:,:,subrun));
        Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
        Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
        [Nc,Nr] = meshgrid(Nc,Nr);
        
        
       Stack_Reg(:,:,subrun) = abs(double(ifft2(fft2(Stack_In(:,:,subrun)).*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc)))));
    
        fprintf('correct cell drift %d of %d',subrun, p )
        fprintf('\n')
    
end

 %Stack_Reg = uint16(Stack_Reg(:,:,:));