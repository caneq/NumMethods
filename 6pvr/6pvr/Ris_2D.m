 Precision = 'double'; 
fidp = fopen ('Param.dat', 'r', 'l'); 
if (fidp == -1) 
   disp('File "Param.dat" not found'); 
    return; 
end 
datap = fread(fidp, 2, 'int'); 
fclose (fidp); 
NX = datap(1); 
NY = datap(2); 
Size = [NX NY]; 
fid = fopen ('Pole.dat', 'r', 'l'); 
if (fid == -1) 
   disp('File "Pole.dat" not found'); 
   return; 
end 
   U = fread (fid, Size, Precision); 
SizeS = size(U); 
x=1:NX; 
y=1:NY; 
[yy, xx] = meshgrid(y,x); 
figure(1);
surf(xx, yy, U) 
xlabel('X') 
ylabel('Y') 
zlabel('U') 
fclose (fid); 