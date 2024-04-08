[fname,pathname]=uigetfile('c:\*.*','Select first file'); %get pathname
filename = strcat(pathname,fname);
image_in=imread(filename,'jpg');
figure;imagesc(image_in);axis equal;title("Original Image of Diwata");
%Now calculate the monochrome luminance by combining the RGB values
%according to the NTSC standard, which applies coefficients related
%to the eye's sensitivity to RGB colors.
gray_im = .2989*image_in(:,:,1)+.5870*image_in(:,:,2)+.1140*image_in(:,:,3);
figure;imagesc(gray_im);colormap(gray);axis equal;title("Original Gray Image");
%disp('grayscale image is stored in gray_im')
%disp('color image is stored in image_in')
my_pic=gray_im; 
[resx,resy] = size(my_pic);
%save the gray_im data to your own variable name
my_picft=fftshift(fft2(my_pic)); 
%calculate the FT of your data
figure;imagesc(log10(1+abs(my_picft)));colormap(gray);title("Fourier Transform of Diwata");
H = fspecial('Gaussian',[resx resy],200);
H=mat2gray(H);
%figure;imagesc(H); colorbar;
my_picft_filt=my_picft.*H; 
figure;imagesc(log10(1+abs(my_picft_filt))); colormap(gray);title("Convoluted Fourier Transform");
my_pic_filt=ifft2(my_picft_filt);
figure;imagesc(abs(my_pic_filt)); colormap(gray);title("Inverse Fourier Transform of Convolution");
H = fspecial('Gaussian',[resx resy],20);
H=mat2gray(H);
%figure;imagesc(H);
my_picft_filt=my_picft.*H;
figure;imagesc(log10(1+abs(my_picft_filt))); colormap(gray);title("Convulted Image"); 
my_pic_filt=ifft2(my_picft_filt); 
%take the inv FT of the new filtered data to create an image
figure;imagesc(abs(my_pic_filt)); colormap(gray); title("Blurred Diwata");