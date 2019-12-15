clc 
clear all 
Precision = 'double'; 
fidp = fopen ('Param.dat', 'r', 'l'); 
if (fidp == -1) 
   disp('File "Param.dat" not found'); 
    return; 
end 
datap = fread (fidp, 3, 'int'); 
fclose (fidp); 
NX = datap(1);
NY = datap(2); 
NT = datap(3); 
Size = [NX NY]; 
fid = fopen ('T1.dat', 'r', 'l'); 
if (fid == -1) 
   disp('File "T1.dat" not found'); 
   return; 
end 
U = fread (fid, Size, Precision); 
SizeS = size(U); 
x=1:NX; 
y=1:NY; 
[yy, xx] = meshgrid(y,x); 
surf(xx, yy, U) 
axis([1 NX 1 NY 0 15]) 
xlabel('X') 
ylabel('Y') 
zlabel('U') 
fclose (fid); 
basename = 'T'; 
for i=2:NT+1 
    filename = sprintf ('%s%d.dat', basename, i); 
    fid = fopen (filename, 'r', 'l'); 
    if (fid == -1) 
        disp('File "T.dat" not found'); 
        return; 
    end 
    U = fread (fid, Size, Precision); 
    SizeS = size(U); 
    x=1:NX; 
    y=1:NY; 
    [yy, xx] = meshgrid(y,x); 
    surf(xx, yy, U) 
    axis([1 NX 1 NY 0 15])
    view(0,45);
    title(['n=',num2str(i),' N=',num2str(NT+1)]) 
    fclose (fid); pause (0.025); 
end