function [h_image,s_image,v_image]= img16rgbtohsv(rgb_image)
    [m,n,c]= size(rgb_image);
    rgb_image=single(rgb_image);
    h_image = zeros(m,n);
    s_image = zeros(m,n);
%     v_image = zeros(m,n);

    mat_maxRGB=single(max(max(rgb_image(:,:,1),rgb_image(:,:,2)),rgb_image(:,:,3)));
    mat_minRGB=single(min(min(rgb_image(:,:,1),rgb_image(:,:,2)),rgb_image(:,:,3)));

    mat_delta=mat_maxRGB-mat_minRGB;

    s_image(mat_maxRGB==0)=0;
    s_image(mat_maxRGB~=0)=mat_delta(mat_maxRGB~=0)./mat_maxRGB(mat_maxRGB~=0);
    
%     if maxRGB == 0
%         hsv(2)=0; %s
%     else
%         hsv(2)=delta/maxRGB;
%     end
    
    t=(rgb_image(:,:,2)-rgb_image(:,:,3))./mat_delta;
    h_image(mat_delta~=0 & mat_maxRGB ==rgb_image(:,:,1))=t(mat_delta~=0 & mat_maxRGB ==rgb_image(:,:,1));

    t=2+(rgb_image(:,:,3)-rgb_image(:,:,1))./mat_delta;
    h_image(mat_delta~=0 & mat_maxRGB ==rgb_image(:,:,2))=t(mat_delta~=0 & mat_maxRGB ==rgb_image(:,:,2));
    
    t=4+(rgb_image(:,:,1)-rgb_image(:,:,2))./mat_delta;
    h_image(mat_delta~=0 & mat_maxRGB ==rgb_image(:,:,3))=t(mat_delta~=0 & mat_maxRGB ==rgb_image(:,:,3));
    
%     if delta~=0
%         if maxRGB == rgb(1) 
%             hsv(1) = (rgb(2)-rgb(3))/delta;
%         elseif maxRGB == rgb(2)
%             hsv(1) = (2+(rgb(3)-rgb(1))/delta);
%         elseif maxRGB == rgb(3)
%             hsv(1) = (4+(rgb(1)-rgb(2))/delta);
%         end       
%     end


    h_image=h_image*60;
    h_image(h_image<0)=h_image(h_image<0)+360;
    
%    
%     hsv(1)=hsv(1)*60;
%     if(hsv(1)<0)
%         hsv(1)=hsv(1)+360;
%     end

    v_image=mat_maxRGB/65535;
    
end