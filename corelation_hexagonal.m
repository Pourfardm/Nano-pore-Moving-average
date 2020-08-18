clear all;
clc;
close all;
% I=imread('C:\Users\office\Pictures\SEM-Arab\1-1-2-.tif');
% I=imread('H:\MATLAB\nano paper\original kafaee paper\SEM-Arab\1-1B-2-.tif');
I=imread('C:\Users\nemat\Pictures\pics\special pics.jpg');
% I=imread('C:\Users\nemat\Pictures\pics\Shingu02-04-315.tif');
I=imresize(I,0.5);
% load I
template_size=17;%17 is the window size of the template
f=@(x,y,xc,yc,sig)(exp(-((x-xc).^2+(y-yc).^2)/sig)); %function handle
x=-template_size:template_size;y=x;
[x,y]=meshgrid(x,y);
threshold=32;%Why 32 ??????????lower bound of correlation
r=5; %by changing r we have different domain size
tic
for theta0=0:59
    s=f(x,y,0,0,5); % Sigma
    for i=1:6
        theta=i*pi/3;
        xc=r*cos(theta+theta0*pi/180);yc=r*sin(theta+theta0*pi/180);
        s=s+f(x,y,xc,yc,5);
    end
%         imshow(uint8(255*s));xlabel('6 lobes Gaussian');  %For rotating the template
%         [XXX, YYY]=meshgrid(1:size(s,1),1:size(s,2));
%         mesh(XXX,YYY,s); colormap('jet'); pause
        
    J(:,:,theta0+1)=xcorr2(I(:,:,1),s)/(template_size^2);
    theta0
end

% figure;
[m_val,k]=max(J,[],3);%maximum on the 3rd dimension
% imshow(uint8(8*min(J,[],3)))
% figure;
% imshow(m>32)
mm=200;
% I=rgb2gray(I);
Kr=mm*trimf(k/60,[0,.25,1])+256-mm;%the same size (x,y axis) of J
Kg=mm*trimf(k/60,[0,.75,1])+256-mm;
Kb=mm*max(trimf(k/60,[-1,0,.5]),...
    trimf(k/60,[.5,1,2]))+256-mm;
I=padarray(I,[(template_size),(template_size)]);%to adapt with Kr
I5=I;
II1=uint8(double(I5(:,:,1)).*Kr/255);
II2=uint8(double(I5(:,:,2)).*Kg/255);
II3=uint8(double(I5(:,:,3)).*Kb/255);
I5(:,:,1)=II1;I5(:,:,2)=II2;I5(:,:,3)=II3;
I5(1:template_size,:,:)=[];I5(end-template_size+1:end,:,:)=[];
I5(:,1:template_size,:)=[];I5(:,end-template_size+1:end,:)=[];
figure;imshow(I5);title(['Without threshold and r=',num2str(r)]);
% labelcolim(:,:,1)=Kr;labelcolim(:,:,2)=Kg;labelcolim(:,:,3)=Kb;
% I4=[];
%  I=padarray(I,[(17),(17)]);
% I4(:,:,1)=I;I4(:,:,2)=I;I4(:,:,3)=I;
% I4=double(I4);
% mrgim=uint8(I4/256.*labelcolim);
I4=I;

II1=I4(:,:,1);
II1(m_val<threshold)=Kr(m_val<threshold);%threshold is lower bound of the corelation
II2=I4(:,:,2);
II2(m_val<threshold)=Kg(m_val<threshold);
II3=I4(:,:,3);
II3(m_val<threshold)=Kb(m_val<threshold);
true_index=m_val<threshold;%index of non-ordered areas
k(true_index)=61;%Assign 61 to non-ordered area

I4(:,:,1)=II1;I4(:,:,2)=II2;I4(:,:,3)=II3;
I4(1:template_size,:,:)=[];I4(end-template_size+1:end,:,:)=[];
I4(:,1:template_size,:)=[];I4(:,end-template_size+1:end,:)=[];
figure;imshow(I4);title('different domains orientation with different colors with threshold');%title('r=9');
%% unsupervised segmentqtion
k(1:template_size,:,:)=[];k(end-template_size+1:end,:,:)=[];%k is index of maximum correlation
k(:,1:template_size,:)=[];k(:,end-template_size+1:end,:)=[];
I(1:template_size,:,:)=[];I(end-template_size+1:end,:,:)=[];
I(:,1:template_size,:)=[];I(:,end-template_size+1:end,:)=[];
for i=0:59
    KK(i+1)=numel(find(k==i));
end
figure,plot(KK);title('Original bin');
bin=KK;
n1=5;
shift1=3;
n2=5;
shift2=3;

filt1=hamming(n1);filt1=filt1/sum(filt1);
bin1=conv(bin,filt1);
bin1=[bin1(shift1:end),bin1(1:shift1-1)];
%     plot(bin1,'r');
% n2=20;


filt2=hamming(n2);filt2=filt2/sum(filt2);
bin2=conv(bin1,filt2);
bin2=[bin2(shift2:end),bin2(1:shift2-1)];
hold on;
plot(bin2,'r');
bias=1;
bin2min=[bin2(end),bin2(1:end-bias)];% one pixel to the right
bin2plus=[bin2(bias+1:end),bin2(bias)];% one pixel to the left
bias=2;
bin2min2=[bin2(end-bias+1:end),bin2(1:end-bias)];% two pixel to the right
bin2plus2=[bin2(bias+1:end),bin2(1:bias)];% two pixel to the left
bias=3;
bin2min3=[bin2(end-bias+1:end),bin2(1:end-bias)];
bin2plus3=[bin2(bias+1:end),bin2(1:bias)];
bias=4;
bin2min4=[bin2(end-bias+1:end),bin2(1:end-bias)];
bin2plus4=[bin2(bias+1:end),bin2(1:bias)];
binfinal=( (bin2>0)...
    &(bin2>=bin2min)   & (bin2>=bin2plus)  & (bin2>=bin2min2) & (bin2>=bin2plus2)...
    & (bin2>=bin2min3) & (bin2>=bin2plus3) & (bin2>=bin2min4) & (bin2>=bin2plus4));
%         & (bin2>=bin2min5) & (bin2>=bin2plus5) & (bin2>=bin2min6) & (bin2>=bin2plus6)...
%         & (bin2>=bin2min7) & (bin2>=bin2plus7) & (bin2>=bin2min8) & (bin2>=bin2plus8));

plot(find(binfinal==1),bin2(binfinal),'mo');
candidate=find(binfinal==1);

%% new label
for iter=1:60
    %     [v INDEX]=min((iter-ind).^2);
    [v INDEX]=min((iter-candidate).^2);
    new_label(iter)=INDEX; %assign each label to its nearest neighbor
end


m=1;
zaviyeh=numel(candidate);
new_label(61)=zaviyeh+1;%non-ordered areas
zaviyeh=zaviyeh+1;

Krr=m*trimf(0:zaviyeh-1,[0,.25,1]*zaviyeh);
Kgg=m*trimf(0:zaviyeh-1,[0,.75,1]*zaviyeh);
Kbb=m*max(trimf(0:zaviyeh-1,[-1,0,.5]*zaviyeh),...
    trimf(0:zaviyeh-1,[.5,1,2]*zaviyeh));% for two part trimf


I6=I;
II1=uint8(double(I6(:,:,1)).*Krr(new_label(k)));
II2=uint8(double(I6(:,:,2)).*Kgg(new_label(k)));
II3=uint8(double(I6(:,:,3)).*Kbb(new_label(k)));
I6(:,:,1)=II1;I6(:,:,2)=II2;I6(:,:,3)=II3;
I6(1:template_size,:,:)=[];I6(end-template_size+1:end,:,:)=[];
I6(:,1:template_size,:)=[];I6(:,end-template_size+1:end,:)=[];
num_domain_reduced=numel(candidate)+1;
figure;imshow(I6);title(['unsupervised segmentation with ',num2str(num_domain_reduced),' domains']);

%% Merge regions of image
k_new=new_label(k);
xt1=num_domain_reduced+1;%number of domains 
xt2=0;%number of domains after merging

domain_num=0;
final_domain=zeros(size(I,1),size(I,2));

for domain=1:xt1
    I_color_test=(k_new==domain);
    % labeling
    [L, num] = bwlabel(I_color_test,4);
    final_domain(L~=0)=final_domain(L~=0)+L(L~=0)+domain_num;
    domain_num=domain_num+num;
%      
end
xt2=domain_num;

m=1;
zaviyeh=domain_num;

Krrr=m*trimf(0:zaviyeh-1,[0,.25,1]*zaviyeh);
Kggg=m*trimf(0:zaviyeh-1,[0,.75,1]*zaviyeh);
Kbbb=m*max(trimf(0:zaviyeh-1,[-1,0,.5]*zaviyeh),...
    trimf(0:zaviyeh-1,[.5,1,2]*zaviyeh));% for two part trimf

I_final(:,:,1)=double(I(:,:,1)).*Krrr(final_domain);
I_final(:,:,2)=double(I(:,:,2)).*Kggg(final_domain);
I_final(:,:,3)=double(I(:,:,3)).*Kbbb(final_domain);
I_final=uint8(I_final);

figure;imshow(I_final),title(['Final domain with ',num2str(domain_num),' domains']);
%% largest domain size
largest_domain_size=0;
for i=1:domain_num
    temp=numel(find(final_domain==i));
if(temp>largest_domain_size)
    largest_domain_size=temp;
    best_i=i;
end
end
largest_domain=(final_domain==best_i);
figure;imshow(largest_domain*255);title(['largest domain size with area of ',num2str(largest_domain_size/numel(final_domain))]);
toc



