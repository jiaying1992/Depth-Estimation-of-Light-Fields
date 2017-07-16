function [mask_ap1,mask_ap2,mask_ap3]=maskap(IM_Pinhole,LF_parameters,dir)
x_size=LF_parameters.x_size;
y_size=LF_parameters.y_size;
UV_diameter=LF_parameters.UV_diameter;
UV_radius=LF_parameters.UV_radius;
UV_size=LF_parameters.UV_size;
mask_ap1=zeros(y_size*UV_diameter,x_size*UV_diameter);
mask_ap2=zeros(y_size*UV_diameter,x_size*UV_diameter);
mask_ap3=zeros(y_size*UV_diameter,x_size*UV_diameter);
% mask_ap4=zeros(y_size*UV_diameter,x_size*UV_diameter);
% mask_ap5=zeros(y_size*UV_diameter,x_size*UV_diameter);
IM_Pinhole1=zeros(y_size+UV_radius*2,x_size+UV_radius*2,3);
IM_Pinhole1(UV_radius+1:y_size+UV_radius,UV_radius+1:x_size+UV_radius,:)=IM_Pinhole;
for y=1:1:y_size
    for x=1:1:x_size
        if(dir(y,x)==0)
            % contrast-peak detection
           mask_ap1(y*UV_diameter-2*UV_radius:y*UV_diameter,x*UV_diameter-2*UV_radius:x*UV_diameter)=1;
           mask_ap2(y*UV_diameter-2*UV_radius:y*UV_diameter,x*UV_diameter-2*UV_radius:x*UV_diameter)=1;
           mask_ap3(y*UV_diameter-2*UV_radius:y*UV_diameter,x*UV_diameter-2*UV_radius:x*UV_diameter)=1;
        else
            mask1=zeros(UV_diameter,UV_diameter);
            mask2=zeros(UV_diameter,UV_diameter);
            mask3=zeros(UV_diameter,UV_diameter);
%             mask4=zeros(UV_diameter,UV_diameter);
%             mask5=zeros(UV_diameter,UV_diameter);
            data=reshape(IM_Pinhole1(y:y+2*UV_radius,x:x+2*UV_radius,:),[UV_size,3]);
            M=UV_size*UV_size-UV_size;
            s=zeros(M,3); %% Make ALL N^2-N similarities
            j=1;
            for i=1:UV_size
                  for k=[1:i-1,i+1:UV_size]
                      s(j,1)=i; 
                      s(j,2)=k; 
                      s(j,3)=-sum((data(i,:)-data(k,:)).^2);
                      j=j+1;
                   end; 
             end;
             p=min(s(:,3)); %% Set preference to median similarity
             [idx,netsim,dpsim,expref]=apcluster(s,p,'plot');
             ide=unique(idx);
             if(length(ide)==1)
                 mask1=1;
                 mask2=1;
                 mask3=1;
             elseif(length(ide)==2)
                 mask1(idx==ide(1))=1;
                 mask2=mask1;
                 mask3=~mask2;
             else 
                 mask1(idx==ide(1))=1;
                 mask2(idx==ide(2))=1;
                 mask3(idx==ide(3))=1;
             end
             centres=idx((UV_size-1)/2+1);
%              mask4(idx==centres)=1;
%              mask5=~mask4;
             mask_ap1(y*UV_diameter-2*UV_radius:y*UV_diameter,x*UV_diameter-2*UV_radius:x*UV_diameter)=mask1;
             mask_ap2(y*UV_diameter-2*UV_radius:y*UV_diameter,x*UV_diameter-2*UV_radius:x*UV_diameter)=mask2;
             mask_ap3(y*UV_diameter-2*UV_radius:y*UV_diameter,x*UV_diameter-2*UV_radius:x*UV_diameter)=mask3;
%              mask_ap4(y*UV_diameter-2*UV_radius:y*UV_diameter,x*UV_diameter-2*UV_radius:x*UV_diameter)=mask4;
%              mask_ap5(y*UV_diameter-2*UV_radius:y*UV_diameter,x*UV_diameter-2*UV_radius:x*UV_diameter)=mask5;
        end
    end
end
end