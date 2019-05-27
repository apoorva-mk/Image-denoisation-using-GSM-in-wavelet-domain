function [yh,yl]=d2wavelet(x,Fs,level)
%y=d2wavelet(x,Fs,level) does the 2nd order Daubchies Wavelet
%Transform of signal x with a sampling frequency 'Fs' and the DWT is
%decomposition is done upto a 'level'
%It returns the matrix of all decompositions and the final approximations.
%
%Instead of using the matlab's inbuilt DWT function, this file explains the
%algorithm for DWT. Mostly useful for learning & academic purposes.
%For other wavelets, the filter values alone can be changed or WFILTERS can
%be used.
%
%The function basically is for Condition Monitoring of rotating equipments
%by vibration based bearing fault diagnosis by the author.
%
%Example:
%       clear all;
%       t=[0:0.0003:8*pi];x=sin(5000*t)+sin(1000*t);
%       x=x(1:2^16);
%       level=5;Fs=1/0.003;
%       d2wavelet(x,Fs,level);
%
% Dont forget to rate or comment on the matlab central site
%http://www.mathworks.in/matlabcentral/fileexchange/authors/258518
%
%https://sites.google.com/site/santhanarajarunachalam/
%Author:Santhana Raj.A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
%Direct Daubchies 2 Wavelet Coeff 
%h1,h2,h3 & h4 are the LPF filter coeffetients for D2
h(1)=(1+sqrt(3))/(4*sqrt(2));
h(4)=(1-sqrt(3))/(4*sqrt(2));
h(2)=(3+sqrt(3))/(4*sqrt(2));
h(3)=(3-sqrt(3))/(4*sqrt(2));
%g1,g2,g3 & g4 are HPF filter coefficients for D4
g(1)=h(4);g(2)=-h(3);g(3)=h(2);g(4)=-h(1);
%Instead of values, one can use the ready matlab function also
%[g,h,g_r,h_r]=wfilters('db2');
x_original=x;
%x=x(1:65536);%algorithm is written for no of values of x, that are multiples of 2
% Making Length of X as a multiple of 2
L=2;level_max=1;
while length(x)>L
L=L*2;
level_max=level_max+1;
end
x=x(1:L);
%estimating level
if level>level_max
    disp('Length of X is small than 2^(level+1). Level taken as max possible value');
    level=level_max-1;
    disp(level);
end
%preallocating for Speed
n=zeros(level+1,1);
yh=zeros(level,length(x)/2);
yl=zeros(level,length(x)/2);
 n(1)=length(x);
for j=1:level
   
    for i=1:n(j)/2-2
        yh(j,i)=x(mod((2*i-1),n(j)))*g(4)+x(mod((2*i),n(j)))*g(3)+x(mod((2*i+1),n(j)))*g(2)+x(mod((2*i+2),n(j)))*g(1);
        yl(j,i)=x(mod((2*i-1),n(j)))*h(4)+x(mod((2*i),n(j)))*h(3)+x(mod((2*i+1),n(j)))*h(2)+x(mod((2*i+2),n(j)))*h(1);
    end
    %for spillover filter variables
    i=n(j)/2-1;
    if(i~=0)
        yh(j,i)=x(mod((2*i-1),n(j)))*g(4)+x(mod((2*i),n(j)))*g(3)+x(mod((2*i+1),n(j)))*g(2)+x((2*i+2))*g(1);
        yl(j,i)=x(mod((2*i-1),n(j)))*h(4)+x(mod((2*i),n(j)))*h(3)+x(mod((2*i+1),n(j)))*h(2)+x((2*i+2))*h(1);
    end
        
    i=n(j)/2;
    yh(j,i)=x(mod((2*i-1),n(j)))*g(4)+x((2*i))*g(3)+x(mod((2*i+1),n(j)))*g(2)+x(mod((2*i+2),n(j)))*g(1);
    yl(j,i)=x(mod((2*i-1),n(j)))*h(4)+x((2*i))*h(3)+x(mod((2*i+1),n(j)))*h(2)+x(mod((2*i+2),n(j)))*h(1);
    
    %reloading x for next loop
    clear x;
    x=yl(j,:);% taking Low pass filter values alone for decomposing to next level
    n(j+1)=n(j)/2;
end
N=2048;T=N/Fs;freq_s=(0:N-1)/T;
%plotting of Original signal
figure();
sig_f=abs(fft(x_original(1:N),N));sig_n=sig_f/(norm(sig_f));plot(freq_s,sig_n);title('Orijinal signal');
%plotting FFTs of all Decompositions
for j=2:level
    clear freq_s;
    N=n(j);
    T=N/Fs;freq_s=(0:N-1)/T;
    figure();
    sig_f=abs(fft(yh(j,1:N),N));sig_n=sig_f/(norm(sig_f));plot(freq_s,sig_n);title(j);
end
%plot FFTs of the last approximation
figure();
sig_f=abs(fft(yl(j,1:N),N));sig_n=sig_f/(norm(sig_f));plot(freq_s,sig_n);title('Approximation');
