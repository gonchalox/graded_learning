% Learn how to get EXR files from TIF (8bits) (tonemapping)
% For each channel
% 1. Create joint histogram (X:TIF, y:EXR)
% 2. Get transformation curve
% 3. Remove outliers and noise
% 4. Fiting to smooth b-spline

% Clear variables
clear;

%%% PARAMETERS %%
% Read from CSV file
exr_dir = '/media/gluzardo/Data/Stuttgart/carousel_fireworks/carousel_fireworks_06/';
tif_dir = '/media/gluzardo/Data/Stuttgart-color-graded/carousel_fireworks/carousel_fireworks_06/';
out_dir = '/media/gluzardo/Data/reg_results_itm/carousel_fireworks/carousel_fireworks_06/';
border=11;
ft=1;
lt=528;
fj=25;
test_list=[1:528];
smooth=0.05;
exp_id='t1';
gamma=1/2.2;
%%%%%%%%%%%%%%%%%

% Add open-exr library 
addpath('/u/gluzardo/Documents/phd/openexr-matlab-master');
exr_files = dir(strcat(exr_dir,'*.exr'));
tif_files = dir(strcat(tif_dir,'*.tif'));      

% Compute minimum and maximum over the whole EXR secuence
mi = 2^16-1;
ma = 0; 

disp('Getting min and max from EXR files..');
for i=ft:lt %All files to normalize% Add open-exr library 
    im = exrread(strcat(exr_dir,exr_files(i).name));
    mi = min(mi,min(im(:)));
    ma = max(ma,max(im(:)));
end
mami= ma - mi;

%Create empty histogram
max_val_y=2^16-1; %max val possible in EXR files
max_val_x=2^8-1; %max val possible in tif files
channels = 3;

r_hist = zeros(max_val_y+1,max_val_x+1);
g_hist = zeros(max_val_y+1,max_val_x+1);
b_hist = zeros(max_val_y+1,max_val_x+1);

disp('Creating joint histogram...');
% To store the maximun value of tif image
ma_tif = 0;
%Joint Hist Calc

for i=ft:fj:lt
      disp(strcat(strcat('Loading TIF:..',tif_files(i).name(1:end-4))));
      disp(strcat(strcat('Loading EXR:.. ',exr_files(i).name(1:end-4))));
      
      exr_image = exrread(strcat(exr_dir,exr_files(i).name));
      tif_image = imread(strcat(tif_dir,tif_files(i).name));

      % Resample to 8 bits
      tif_image=uint8((2^8-1)*(single(tif_image)/(2^16-1)));
      
      %Max tif value in samples, it used to show the curve
      ma_tif = max(ma_tif,max(tif_image(:)));
            
      exr_data = (exr_image - mi)./mami;
      exr_data = uint16((2^16-1)*(exr_data.^gamma));

      %Accum in joint hist per channel
      for m=border+1:size(exr_data,1)-border
        for n=border+1:size(exr_data,2)-border
            y=exr_data(m,n,1)+1;
            x=tif_image(m,n,1)+1;
            r_hist(y,x)=r_hist(y,x)+1;
            
            y=exr_data(m,n,2)+1;
            x=tif_image(m,n,2)+1;
            g_hist(y,x)=g_hist(y,x)+1;
             
            y=exr_data(m,n,3)+1;
            x=tif_image(m,n,3)+1;
            b_hist(y,x)=b_hist(y,x)+1;
        end
      end  
end

%Visualize
ji_hr=figure;
image(r_hist);
xlabel('TIF','fontsize',10);
ylabel('EXR','fontsize',10);
title('Red channel','fontsize',10);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
saveas(ji_hr,strcat(out_dir,strcat(exp_id,'_r_joint_hist.fig')),'fig')

ji_hg=figure;
image(g_hist);
xlabel('TIF','fontsize',10);
ylabel('EXR','fontsize',10);
title('Green channel','fontsize',10);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
saveas(ji_hg,strcat(out_dir,strcat(exp_id,'_r_joint_hist.fig')),'fig')

ji_hb=figure;
image(b_hist);
xlabel('TIF','fontsize',10);
ylabel('EXR','fontsize',10);
title('Blue channel','fontsize',10);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
saveas(ji_hb,strcat(out_dir,strcat(exp_id,'_r_joint_hist.fig')),'fig')


%Save files
%disp('Save results to file...');
%save(strcat(out_dir,'results_full'),'r_hist','g_hist','b_hist');
% save('results_bins','r_bin','g_bin','b_bin');

%Row data curves
disp('Creating row data dataset...');
n_data = length(r_hist(r_hist>0));
xr=zeros(n_data,1);
yr=zeros(n_data,1);
xg=zeros(n_data,1);
yg=zeros(n_data,1);
xb=zeros(n_data,1);
yb=zeros(n_data,1);

%Calcule median datasets  R
disp('Creating median dataset..');
median_rx=zeros(1,2^8);
median_ry=zeros(1,2^16);
cont=1;
for i=1:2^8
    s = find(r_hist(:,i)>0);
    if(~isempty(s))
        median_ry(cont)=median(s);
        median_rx(cont)=i;
        cont=cont+1;
    end
end
median_rx = median_rx(1:cont-1);
median_ry = median_ry(1:cont-1);

%Calcule median G
median_gx=zeros(1,2^8);
median_gy=zeros(1,2^16);
cont=1;
for i=1:2^8
    s = find(g_hist(:,i)>0);
    if(~isempty(s))
        median_gy(cont)=median(s);
        median_gx(cont)=i;
        cont=cont+1;
    end
end
median_gx = median_gx(1:cont-1);
median_gy = median_gy(1:cont-1);

%Calcule median B 
median_bx=zeros(1,2^8);
median_by=zeros(1,2^16);
cont=1;
for i=1:2^8
    s = find(b_hist(:,i)>0);
    if(~isempty(s))
        median_by(cont)=median(s);
        median_bx(cont)=i;
        cont=cont+1;
    end
end
median_bx = median_bx(1:cont-1);
median_by = median_by(1:cont-1);
% 
% %Clear mem
% % clear('r_hist');
% % clear('g_hist');
% % clear('b_hist');
% 
% %Save datasets
% %dataset = [median_rx median_ry];
% % dlmwrite('median_datasetr.dat',dataset,'delimiter',',');
% % dataset = [xg yg];
% % dlmwrite('median_datasetg.dat',dataset,'delimiter',',');
% % dataset = [xb yb];
% % dlmwrite('median_datasetb.dat',dataset,'delimiter',',');
% 
% 

%Remove outliers
disp('Remove outliers...');
tol=8000;
median_ry=[0,median_ry];
median_rx=[-1,median_rx];
ind=median_rx;
for i=2:length(median_ry)
    if(abs(median_ry(i-1)-median_ry(i))>tol)
        ind(i)=-1;
        median_ry(i)=median_ry(i-1);
    end    
end
median_rx=(median_rx(ind>=0));
median_ry=(median_ry(ind>=0));

median_gy=[0,median_gy];
median_gx=[-1,median_gx];
ind=median_gx;
for i=2:length(median_gy)
    if(abs(median_gy(i-1)-median_gy(i))>tol)
        ind(i)=-1;
        median_gy(i)=median_gy(i-1);
    end    
end
median_gx=(median_gx(ind>=0));
median_gy=(median_gy(ind>=0));

median_by=[0,median_by];
median_bx=[-1,median_bx];
ind=median_bx;
for i=2:length(median_by)
    if(abs(median_by(i-1)-median_by(i))>tol)
        ind(i)=-1;
        median_by(i)=median_by(i-1);
    end    
end
median_bx=(median_bx(ind>=0));
median_by=(median_by(ind>=0));

%Remove noise
disp('Remove noise...');
window_size=5;
k=ones(1,window_size)/window_size;
median_ry = filter(k,1,median_ry);
median_gy = filter(k,1,median_gy);
median_by = filter(k,1,median_by);

% Add start and end point
median_ry=[0,median_ry,2^16-1];
median_rx=[0,median_rx,2^8-1];

median_gy=[0,median_gy,2^16-1];
median_gx=[0,median_gx,2^8-1];

median_by=[0,median_by,2^16-1];
median_bx=[0,median_bx,2^8-1];

%Modify curve before fitting
disp('Adjusting curve...');
median_ry=[-median_ry(end:-1:1), median_ry];
sp=median_rx-median_rx(1);
sp = -sp(end:-1:1)+median_rx(1);
median_rx=[sp,median_rx];

median_gy=[-median_gy(end:-1:1), median_gy];
sp=median_gx-median_gx(1);
sp = -sp(end:-1:1)+median_gx(1);
median_gx=[sp,median_gx];

median_by=[-median_by(end:-1:1), median_by];
sp=median_bx-median_bx(1);
sp = -sp(end:-1:1)+median_bx(1);
median_bx=[sp,median_bx];

% Fitting
disp('Fitting smooth spline...');
ppr=csaps(median_rx,median_ry,smooth);
hr=figure;
plot(median_rx,median_ry,'ok'); hold on; fnplt(ppr,[0,2^8]); hold off
axis([0,2^8,0,2^16])
saveas(hr,strcat(out_dir,strcat(exp_id,'_r.fig')),'fig')
   
ppg=csaps(median_gx,median_gy,smooth);
hg=figure;
plot(median_gx,median_gy,'ok'); hold on; fnplt(ppg,[0,2^8]); hold off
axis([0,2^8,0,2^16])
saveas(hg,strcat(out_dir,strcat(exp_id,'_g.fig')),'fig')

ppb=csaps(median_bx,median_by,smooth);
hb=figure;
plot(median_bx,median_by,'ok'); hold on; fnplt(ppb,[0,2^8]); hold off
axis([0,2^8,0,2^16])
saveas(hb,strcat(out_dir,strcat(exp_id,'_b.fig')),'fig')

% %Creating lookup table
% lu_table=zeros(2^16,4);
% lu_table(:,1)=0:2^16-1;
% lu_table(:,2)=uint16(fnval(ppr,lu_table(:,1)));
% lu_table(:,3)=uint16(fnval(ppg,lu_table(:,1)));
% lu_table(:,4)=uint16(fnval(ppb,lu_table(:,1)));
% dlmwrite(strcat(out_dir,strcat(exp_id,'_loopuptable.dat')),lu_table,'delimiter',',');

disp('Creating automatic color graded images...');

for i=test_list
   % Read EXR and generate a color graded TIF and EXR
    tif_image = imread(strcat(tif_dir,tif_files(i).name));
    tif_image = single((2^8-1)*single(tif_image)/(2^16-1));
    exr_data_itm = zeros(size(tif_image));

    disp(strcat('Grading TIF:..',exr_files(i).name));  
    exr_data_itm(border+1:end-border,border+1:end-border,1)=uint16(fnval(ppr,tif_image(border+1:end-border,border+1:end-border,1)));
    exr_data_itm(border+1:end-border,border+1:end-border,2)=uint16(fnval(ppg,tif_image(border+1:end-border,border+1:end-border,2)));
    exr_data_itm(border+1:end-border,border+1:end-border,3)=uint16(fnval(ppb,tif_image(border+1:end-border,border+1:end-border,3)));
    
    exr_image=((single(exr_data_itm)/(2^16-1)).^(1/gamma))*ma-mi;
    exrwrite(exr_image,strcat(strcat(out_dir,'exr_itm/itm_exr_'),exr_files(i).name));
end





