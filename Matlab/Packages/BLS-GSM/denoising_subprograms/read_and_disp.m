function im = read_and_disp(path)

im = imread('barco.png');
im = double(im);
figure(1)
%rang0 = showIm(im,'auto');title('Original image');