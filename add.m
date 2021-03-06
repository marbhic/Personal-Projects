                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %jpgread2gray is below with whatever image stored in img: 
[fname,pathname]=uigetfile('c:\*.*','Select first file');  %get pathname

filename = strcat(pathname,fname);

image_in=imread(filename,'jpg');

figure;imagesc(image_in);axis equal;title("Original Image of Sparks");
%Now calculate the monochrome luminance by combining the RGB values 
%according to the NTSC standard, which applies coefficients related 
%to the eye's sensitivity to RGB colors. 
gray_im = .2989*image_in(:,:,1)+.5870*image_in(:,:,2)+.1140*image_in(:,:,3);
img = gray_im;
[fname,pathname]=uigetfile('c:\*.*','Select first file');  %get pathname

filename = strcat(pathname,fname);

image_in=imread(filename,'jpg');

figure;imagesc(image_in);axis equal;title("Original Image of Sparks");
%Now calculate the monochrome luminance by combining the RGB values 
%according to the NTSC standard, which applies coefficients related 
%to the eye's sensitivity to RGB colors. 
gray_im3 = .2989*image_in(:,:,1)+.5870*image_in(:,:,2)+.1140*image_in(:,:,3);
img3 = gray_im3;

figure;imagesc(gray_im3);colormap(gray);axis equal;title("Original Gray Image");

[fname,pathname]=uigetfile('c:\*.*','Select first file');  %get pathname

filename = strcat(pathname,fname);

image2=imread(filename,'jpg');

figure;imagesc(image2);axis equal;title("Original Image of Sparks");
%Now calculate the monochrome luminance by combining the RGB values 
%according to the NTSC standard, which applies coefficients related 
%to the eye's sensitivity to RGB colors. 
gray_im2 = .2989*image2(:,:,1)+.5870*image2(:,:,2)+.1140*image2(:,:,3);
%disp('grayscale image is stored in gray_im')
%disp('color image is stored in image_in')


% begin linearity code 

img2=gray_im2;
imgft2=fftshift(fft2(img2));

%ft of image
figure;imagesc(log10(1+abs(imgft2)));colormap(gray);title("Fourier Transform")

[resx1,resy1]=size(img2);
[resx3,resy3]=size(img3);
disp([resx1,resy1])
disp([resx,resy])


addition=double(mat2gray(img))+double(mat2gray(img3))+ double (mat2gray(img2));        
%add the image and the spatial frequency - notes 1) neeeded to convert the smiley to double in this case. 2) needed to decrease dimensions by one on the frequency 3) **need to be sure everything is scaled the same (mat2gray)
              
                                                                  
figure;imagesc(addition);colormap(gray); title("Addition of Images");%display addition image

additionft=fftshift(fft2(addition));  

figure;imagesc(log10(1+abs(additionft)));colormap(gray); title("Fourier Transform of Addition") %FT of addition image