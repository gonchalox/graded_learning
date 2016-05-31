% Clear cariables
clear;

%%% PARAMETERS %%
% Read from CSV file
exr_dir = '/media/gluzardo/Data/Stuttgart/showgirl_02/';
tif_dir = '/media/gluzardo/Data/Stuttgart-color-graded/showgirl_02/';
out_dir = '/media/gluzardo/Data/reg_results/showgirl_02/';
border=15;
ft=1;
lt=156; %300
fj=12;
test_list=[100];
smooth=1e-10;
exp_id='t1';
%%%%%%%%%%%%%%%%%

% Add open-exr library 
addpath('/u/gluzardo/Documents/phd/openexr-matlab-master');


exr_files = dir(strcat(exr_dir,'*.exr'));
tif_files = dir(strcat(tif_dir,'*.tif'));      

% Compute minimum and maximum over the whole EXR secuence
mi = 2^16;
ma = 0; 

disp('Phase 1: Get min and max from EXR files ');
for i=ft:lt %All files to normalize% Add open-exr library 
    im = exrread(strcat(exr_dir,exr_files(i).name));
    mi = min(mi,min(im(:)));
    ma = max(ma,max(im(:)));
end
mami= ma - mi;

%Create empty histogram
min_val=uint16(0);
max_val=uint16(2^16-1);
channels = 3;
r_hist = zeros(2^16,2^16);
g_hist = zeros(2^16,2^16);
b_hist = zeros(2^16,2^16);

% To store the maximun value of tif image
ma_tif = 0;
%Joint Hist Calc
gamma=1/2.2;
for i=ft:fj:lt
      disp(strcat(strcat('Phase 2: Loading EXR:.. ',exr_files(i).name(1:end-4))));
      disp(strcat(strcat('Phase 2: Loading TIF:..',tif_files(i).name(1:end-4))));
      exr_image = exrread(strcat(exr_dir,exr_files(i).name));
      tif_image = imread(strcat(tif_dir,tif_files(i).name));
      
      %Max tif value in samples, it used to show the curve
      ma_tif = max(ma_tif,max(tif_image(:)));
      
      exr_data = (exr_image - mi)./mami;
      exr_data = uint16((2^16-1)*(exr_data.^gamma));

      %Accum in joint hist per channel
      for m=border+1:size(exr_data,1)-border
        for n=border+1:size(exr_data,2)-border
            x=exr_data(m,n,1)+1;
            y=tif_image(m,n,1)+1;
            r_hist(y,x)=r_hist(y,x)+1;
            
            x=exr_data(m,n,2)+1;
            y=tif_image(m,n,2)+1;
            g_hist(y,x)=g_hist(y,x)+1;
             
            x=exr_data(m,n,3)+1;
            y=tif_image(m,n,3)+1;
            b_hist(y,x)=b_hist(y,x)+1;
        end
      end  
end

% % Bins for visualization only
% 
% disp('Creating joint histogram for visualization...');
% bin_size = 16;
% r_bin = zeros((max_val+1)/bin_size,(max_val+1)/bin_size);
% g_bin = zeros((max_val+1)/bin_size,(max_val+1)/bin_size);
% b_bin = zeros((max_val+1)/bin_size,(max_val+1)/bin_size);
% 
% for m=1:(max_val+1)/bin_size-1
%     for n=1:max_val/bin_size
%         reg = r_hist((m-1)*bin_size+1:m*bin_size,(n-1)*bin_size+1:n*bin_size);
%         geg = g_hist((m-1)*bin_size+1:m*bin_size,(n-1)*bin_size+1:n*bin_size);
%         beg = b_hist((m-1)*bin_size+1:m*bin_size,(n-1)*bin_size+1:n*bin_size);
%         
%         r_bin(m,n)=sum(sum(reg)); 
%         g_bin(m,n)=sum(sum(geg)); 
%         b_bin(m,n)=sum(sum(beg)); 
%     end
% end
% 
% figure;
% image(r_bin);axis([0 max_val/bin_size 0 floor(ma_tif/bin_size)]);
% ax=gca;
% xlabel('EXR','fontsize',10);
% ylabel('Gradded TIF','fontsize',10);
% title('Red channel','fontsize',10);
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
% 
% figure;
% image(g_bin);axis([0 max_val/bin_size 0 floor(ma_tif/bin_size)]);
% xlabel('EXR','fontsize',10);
% ylabel('Gradded TIF','fontsize',10);
% title('Green channel','fontsize',10);
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
% 
% figure;
% image(b_bin);axis([0 max_val/bin_size 0 floor(ma_tif/bin_size)]);
% xlabel('EXR','fontsize',10);
% ylabel('Gradded TIF','fontsize',10);
% title('Blue channel','fontsize',10);
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);

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

cont=1;
%Red channel
for n=1:2^16
   for m=1:2^16
       if(r_hist(m,n)>0)
           xr(cont)=m-1;
           yr(cont)=n-1;
           cont=cont+1;
       end
   end 
end
%Green channel
cont=1;
for n=1:2^16
   for m=1:2^16
       if(g_hist(m,n)>0)
           xg(cont)=m-1;
           yg(cont)=n-1;
           cont=cont+1;
       end
   end 
end
%Blue channel
cont=1;
for n=1:2^16
   for m=1:2^16
       if(b_hist(m,n)>0)
           xb(cont)=m-1;
           yb(cont)=n-1;
           cont=cont+1;
       end
   end 
end

%Save datasets
% dataset = [xr yr];
% dlmwrite(strcat(out_dir,'datasetr.dat'),dataset,'delimiter',',');
% dataset = [xg yg];
% dlmwrite(strcat(out_dir,'datasetg.dat'),dataset,'delimiter',',');
% dataset = [xb yb];
% dlmwrite(strcat(out_dir,'datasetb.dat'),dataset,'delimiter',',');


%Calcule median datasets  R
disp('Creating median dataset..');
median_rx=zeros(1,2^16);
median_ry=zeros(1,2^16);
cont=1;
for i=1:2^16
    s = find(r_hist(:,i)>0);
    if(~isempty(s)>0)
        median_ry(cont)=median(s);
        median_rx(cont)=i;
        cont=cont+1;
    end
end
median_rx = median_rx(1:cont-1);
median_ry = median_ry(1:cont-1);

%Calcule median G
median_gx=zeros(1,2^16);
median_gy=zeros(1,2^16);
cont=1;
for i=1:2^16
    s = find(g_hist(:,i)>0);
    if(~isempty(s)>0)
        median_gy(cont)=median(s);
        median_gx(cont)=i;
        cont=cont+1;
    end
end
median_gx = median_gx(1:cont-1);
median_gy = median_gy(1:cont-1);

%Calcule median B 
median_bx=zeros(1,2^16);
median_by=zeros(1,2^16);
cont=1;
for i=1:2^16
    s = find(b_hist(:,i)>0);
    if(~isempty(s)>0)
        median_by(cont)=median(s);
        median_bx(cont)=i;
        cont=cont+1;
    end
end
median_bx = median_bx(1:cont-1);
median_by = median_by(1:cont-1);

%Save datasets
%dataset = [median_rx median_ry];
% dlmwrite('median_datasetr.dat',dataset,'delimiter',',');
% dataset = [xg yg];
% dlmwrite('median_datasetg.dat',dataset,'delimiter',',');
% dataset = [xb yb];
% dlmwrite('median_datasetb.dat',dataset,'delimiter',',');

%  Cubic spline fit

disp('Performing cube spline fit...');
%SMOOTH
ppr=csaps(median_rx,median_ry,smooth);
hr=figure;
plot(median_rx,median_ry,'ok'); hold on; fnplt(ppr); hold off
saveas(hr,strcat(out_dir,strcat(exp_id,'_r.fig')),'fig')
   
ppg=csaps(median_gx,median_gy,smooth);
hg=figure;
plot(median_gx,median_gy,'ok'); hold on; fnplt(ppg); hold off
saveas(hg,strcat(out_dir,strcat(exp_id,'_g.fig')),'fig')

ppb=csaps(median_bx,median_by,smooth);
hb=figure;
plot(median_bx,median_by,'ok'); hold on; fnplt(ppb); hold off
saveas(hb,strcat(out_dir,strcat(exp_id,'_b.fig')),'fig')

%%TOLERANCE
% tol=100000000000;
% ppr=spaps(median_rx,median_ry,tol);
% figure;
% plot(median_rx,median_ry,'ok'); hold on; fnplt(ppr); hold off
%    
% ppg=spaps(median_gx,median_gy,tol);
% figure;
% plot(median_gx,median_gy,'ok'); hold on; fnplt(ppg); hold off
% 
% ppb=spaps(median_bx,median_by,tol);
% figure;
% plot(median_bx,median_by,'ok'); hold on; fnplt(ppb); hold off


%Creating lookup table
lu_table=zeros(2^16,4);
lu_table(:,1)=0:2^16-1;
lu_table(:,2)=uint16(fnval(ppr,lu_table(:,1)));
lu_table(:,3)=uint16(fnval(ppg,lu_table(:,1)));
lu_table(:,4)=uint16(fnval(ppb,lu_table(:,1)));
dlmwrite(strcat(out_dir,strcat(exp_id,'_loopuptable.dat')),lu_table,'delimiter',',');

disp('Creating automatic color graded images...');

%for i=1:nr
for i=test_list
    %Evaluate in a test image
    %Read EXR and generate a color graded TIF and EXR
    exr_image = exrread(strcat(exr_dir,exr_files(i).name));
    exr_data = (exr_image - mi)./mami; %Adjust
    exr_data = ((2^16-1)*(exr_data.^gamma)); %Gamma correction
    tif_image_graded_data = uint16(zeros(size(exr_data)));

    disp(strcat('Grading EXR:..',exr_files(i).name));
    tif_image_graded_data(border+1:end-border,border+1:end-border,1)=uint16(fnval(ppr,exr_data(border+1:end-border,border+1:end-border,1)));
    tif_image_graded_data(border+1:end-border,border+1:end-border,2)=uint16(fnval(ppg,exr_data(border+1:end-border,border+1:end-border,2)));
    tif_image_graded_data(border+1:end-border,border+1:end-border,3)=uint16(fnval(ppb,exr_data(border+1:end-border,border+1:end-border,3)));
    
    % Save the color gradded
    imwrite(tif_image_graded_data,strcat(strcat(out_dir,'automatic_graded_tif/automatic_graded_'),tif_files(i).name));
    % Save original ungraded
    imwrite(uint16(exr_data),strcat(strcat(out_dir,'automatic_graded_tif/ungraded_'),tif_files(i).name));
    % Save theoriginal graded
    tif_image=imread(strcat(tif_dir,tif_files(i).name));
    imwrite(tif_image,strcat(strcat(out_dir,'automatic_graded_tif/graded_'),tif_files(i).name));
    
    %Convert to EXR format and save
    exr_image_graded=((single(tif_image_graded_data)/2^16).^1/gamma)*mami+mi;
    % Save the color graded
    exrwrite(exr_image_graded,strcat(strcat(out_dir,'automatic_graded_exr/automatic_graded_'),exr_files(i).name)); 
    % Save original ungraded
    exrwrite(exr_image,strcat(strcat(out_dir,'automatic_graded_exr/ungraded_'),exr_files(i).name)); 
end




