I = im2double(imread('/eno/cllee3/DATA/esbilili/20150731reprocessed/img0011adj.jpg'));
y = mean(I, 2);
imshow(y)
sigma = 70; % choosen by visual inspection
G = fspecial('gaussian', 3*sigma+1, sigma);
yb = imfilter(I, G, 'replicate');
imshow(yb)
Ic = bsxfun(@minus, I, yb);
subplot(1,2,1), imshow(I)
subplot(1,2,2), imshow(Ic)