function [hsv]= rgb16tohsv16(rgb)
    hsv=[0,0,0];
    maxRGB=single(max(rgb));
    minRGB=single(min(rgb));
    delta=maxRGB-minRGB;
    if maxRGB == 0
        hsv(2)=0; %s
    else
        hsv(2)=delta/maxRGB;
    end    
    if delta~=0
        if maxRGB == rgb(1) 
            hsv(1) = (rgb(2)-rgb(3))/delta;
        elseif maxRGB == rgb(2)
            hsv(1) = (2+(rgb(3)-rgb(1))/delta);
        elseif maxRGB == rgb(3)
            hsv(1) = (4+(rgb(1)-rgb(2))/delta);
        end       
    end
    hsv(1)=hsv(1)*60;
    if(hsv(1)<0)
        hsv(1)=hsv(1)+360;
    end
    hsv(3)=maxRGB/65535;
end