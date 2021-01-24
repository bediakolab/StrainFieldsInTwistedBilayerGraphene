function sigRecon = TGV(sigNoisy,numIter,padNum,alpha,beta,theta1,theta2)

% This Demo repeats the results of the following papers
% Please cite them if you find them useful for your research

% Reference:
% 1: W Lu, J Duan, Z Qiu, Z Pan, W Ryan Liu, L Bai
%    Implementation of high order variational models made easy for image processing

% 2: J Duan, Z Qiu, W Lu, G Wang, Z Pan, L Bai
%    An edge-weighted second order variational model for image decomposition

% code Writen by
% Wenqi Lu and Jinming Duan
% contact email: j.duan@imperial.ac.uk
% March 2018

% Retrieved from https://github.com/j-duan/HOVM,
% Slightly modified by Colin Ophus - 2020 June
% This implementation is based on from ref 1

flagPlot = false;
% if nargin < 3
%     padNum=10;
%     alpha=30;
%     beta=3*alpha;
%     theta1=5;
%     theta2=theta1;
% end

f0=padarray(sigNoisy,[padNum,padNum],'symmetric');
f0=double(f0);
[m,n]=size(f0);

p1=zeros(m,n);
p2=zeros(m,n);
w1=zeros(m,n);
w2=zeros(m,n);
b1=zeros(m,n);
b2=zeros(m,n);

v11=zeros(m,n);
v22=zeros(m,n);
v3 =zeros(m,n);

d11=zeros(m,n);
d22=zeros(m,n);
d3 =zeros(m,n);


[Y,X]=meshgrid(0:n-1,0:m-1);
G1=cos(2*pi*Y/n)-1;
G2=cos(2*pi*X/m)-1;
a11=theta1-2*theta2*(G1+0.5*G2);
a22=theta1-2*theta2*(0.5*G1+G2);
i=sqrt(-1);
a12=-0.5*theta2*(-1+cos(2*pi*X/m)+i*sin(2*pi*X/m))...
    .*(1-cos(2*pi*Y/n)+i*sin(2*pi*Y/n));
a21=-0.5*theta2*(-1+cos(2*pi*Y/n)+i*sin(2*pi*Y/n))...
    .*(1-cos(2*pi*X/m)+i*sin(2*pi*X/m));
D=(theta1-2*theta2*(G1+G2)).*(theta1-theta2*(G1+G2));

if flagPlot == true
    sigMean = median(sigNoisy(:));
    sigRMS = sqrt(mean((sigNoisy(:) - sigMean(:)).^2));
end

for step = 1:numIter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update u using FFT
    div_wpb=Bx(w1+p1-b1)+By(w2+p2-b2);
    g=f0-theta1*div_wpb;
    u=ifft2(fft2(g)./(1-2*theta1*(G1+G2)),'symmetric');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h1=theta1*(Fx(u)+b1-w1)-theta2*(Fx(v11-d11)+Fy(v3-d3));
    h2=theta1*(Fy(u)+b2-w2)-theta2*(Fy(v22-d22)+Fx(v3-d3));
    Fh1=fft2(h1);
    Fh2=fft2(h2);
    p1=ifft2((a22.*Fh1-a12.*Fh2)./D,'symmetric');
    p2=ifft2((a11.*Fh2-a21.*Fh1)./D,'symmetric');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C1=Fx(u)-p1+b1;
    C2=Fy(u)-p2+b2;
    abs_C=sqrt(C1.^2+C2.^2+eps);
    w1=max(abs_C-alpha/theta1,0).*C1./abs_C;
    w2=max(abs_C-alpha/theta1,0).*C2./abs_C;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update w using soft thresholding
    S1=Bx(p1)+d11;
    S2=By(p2)+d22;
    S3=0.5*(By(p1)+Bx(p2))+d3;
    abs_S=sqrt(S1.^2+S2.^2+2*S3.^2+eps);
    v11=max(abs_S-beta/theta2,0).*S1./abs_S;
    v22=max(abs_S-beta/theta2,0).*S2./abs_S;
    v3 =max(abs_S-beta/theta2,0).*S3./abs_S;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update Bregman iterative parameters b
    b1 =C1-w1;
    b2 =C2-w2;
    d11=S1-v11;
    d22=S2-v22;
    d3 =S3-v3;
    
    % plotting
    if flagPlot == true
        if mod(log(step)/log(2),1) == 0 || step == numIter
            figure(11)
            clf
            imagesc([sigNoisy u(padNum+1:m-padNum,padNum+1:n-padNum)])
            text(10,20,[num2str(step) ' / ' num2str(numIter)],...
                'fontsize',16,'color','k')
            axis equal off
            colormap(jet(256))
            caxis([-3 3]*sigRMS + sigMean)
            set(gca,'position',[0 0 1 1])
            drawnow;
        end
    end
end
% toc
sigRecon=u(padNum+1:m-padNum,padNum+1:n-padNum);
% figure; imagesc(u); colormap(gray); axis off; axis equal;
end

% Forward derivative operator on x with periodic boundary condition
function Fxu = Fx(u)
Fxu = circshift(u,[0 -1])-u;
end

% Forward derivative operator on y with periodic boundary condition
function Fyu = Fy(u)
Fyu = circshift(u,[-1 0])-u;
end

% Backward derivative operator on x with periodic boundary condition
function Bxu = Bx(u)
Bxu = u - circshift(u,[0 1]);
end

% Backward derivative operator on y with periodic boundary condition
function Byu = By(u)
Byu = u - circshift(u,[1 0]);
end