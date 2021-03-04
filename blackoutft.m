[fname,pathname]=uigetfile('c:\*.*','Select first file'); %get pathname
filename = strcat(pathname,fname);
image_in=imread(filename,'jpg');
figure;imagesc(image_in);axis equal;title("Original Image of Luyong");
%Now calculate the monochrome luminance by combining the RGB values
%according to the NTSC standard, which applies coefficients related
%to the eye's sensitivity to RGB colors.
gray_im = .2989*image_in(:,:,1)+.5870*image_in(:,:,2)+.1140*image_in(:,:,3);
figure;imagesc(gray_im);colormap(gray);axis equal;title("Original Gray Image");
%disp('grayscale image is stored in gray_im');
%disp('color image is stored in image_in');

beach=gray_im; %creates matrix “beach” from gray_im data
beachft=fft2(beach); %takes 2D FT of beach.
beachft=fftshift(beachft); %brings the zero-frequency data to the origin
beachftlog=log10(1+abs(beachft)); %scales the FT so that you can see all the data and not just a bright spot in the middle
figure;imagesc(beachftlog);colormap(gray); title("Fourier Transform of Luyong"); %create a figure, show the log version of the data, change the colormap to gray
beachftbo=blackout(beachft); %***CALLS THE FUNCTION “blackoutft”*** needs FT data input. Stores output in new variable beachftbo
newbeach=ifft2(beachftbo); %take the inverse FT of the blacked out FT data from the function
figure;imagesc(abs(newbeach));colormap(gray);axis equal; title ("Inverse FT of Blackout Image");


function [imnew] = blackout(imin)
figure;imagesc(log10(1+abs(imin)));colormap(gray); %just displaying
[m,n]=ginput(2); %take in top-left and bottom-right of blackout region from user
m=round(m);
n=round(n);
imnew=imin;
imnew(n(1):n(2),m(1):m(2))=0.0;
figure;imagesc(log10(1+abs(imnew)));colormap(gray);title("Blackout Image of Luyong"); %just displaying
end
