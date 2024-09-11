%%Code to generate the wave heights
%%Spectral analysis of waterlevel inputs
close all
clear all

ncfile='C:\Users\eday\Documents\XBeach\NWFW_NC_Oriental\Trial_Models\trial_model9\xboutput.nc';
%ncinfo=ncinfo(ncfile);
t=ncread(ncfile,'globaltime');
x=ncread(ncfile,'globalx');
zs=ncread(ncfile,'zs');
zs=squeeze(zs)';
zb=ncread(ncfile,'zb');
zb=squeeze(zb)';
x=x',

%wave record 15 minute time period is 900sec.
delta_t=floor(length(zs)/900)

ts=1
for i=1:delta_t
   % n1(i)=ts;
    n=ts;
    ts=ts+899;
    zs_i=zs(n:ts,:);
    zs_int{i}=zs_i;
    zb_i=zb(n:ts,:);
    zb_int{i}=zb_i;
    t_i=t(n:ts);
    t_int{i}=t_i;
end

index_j=length(zs_int);
for j=1:index_j
    zs_i2=zs_int{j};
    zb_i2=zb_int{j};
    zs_i3=zs_i2;
    zs_i2(isnan(zs_i2)) = zb_i2(isnan(zs_i2));
    t_i2=t_int{j};
    
for k = 1:size(zs_i2,2)% for every location
            [f,Snn,A,ff,Af] = comp_spec(detrend(zs_i2(:,k)),1/(t_i2(2)-t_i2(1)));
            df = f(2)-f(1);

           % calculate variance below (IG) and above (short waves) fsplit
            fsplit = 0.05;
            m0hf(k) = sum(Snn(f >= fsplit)*df);% short wave energy density variance
            m0lf(k) = sum(Snn(f <  fsplit)*df);%  long wave energy density variance
                
            Hm0_hfi = 4*sqrt(m0hf);
            Hm0_lfi = 4*sqrt(m0lf);
            Hm0_i = sqrt(Hm0_hfi.^2+Hm0_lfi.^2);
            zs_regi(k)=zs_i3(1,k);
            zs_avgi(k)=mean(zs_i3(:,k));
            zs_avgi2=zs_avgi;
end
Hm0_hf{j}=Hm0_hfi;
Hm0_lf{j}=Hm0_lfi;
Hm0{j}=Hm0_i;

zs_avg{j}=zs_avgi2;
zs_reg{j}=zs_regi;
end
bed=zb(1,1:end);
t_int=[];
for jj=1:delta_t
    n=jj
    t1=n*900
    tstr=datestr(seconds(t1),'DD HH:MM:SS')
    t_int{jj}=tstr
end
save('Transect1_trial9.mat','bed','x','Hm0','zs_avg','zs_reg','t_int','zs');
%%
close all
clear all
[num dat all]=xlsread('jonswap_100yr.xls');
clear dat all
wh_js=num(1:380,1)';

load Transect1_trial9.mat
figure(1)
axis tight manual
set(gca)
%set(gcf, 'Position', get(0, 'Screensize'));
for ii=1:length(Hm0)
    t=t_int{1,ii}
    p1=plot(x,bed,'k')
    hold on
    p2=plot(x,Hm0{1,ii},'r')
    hold on
    p3=plot(x,zs_avg{1,ii},'g')
    hold on 
    p4=plot(x,zs_reg{1,ii},'b')
    hold on 
    p5=plot(-1.345073035543540e+02,wh_js(ii),'ro')
    title_n=strcat('Time',{' '},t);
    title(title_n)
    ylim([-3 3])
    xlim([0 250])
    legend([p1 p2 p3 p4 p5],{'Bedlevel','Hm0','Mean Water Level','Instaneous Water Level','Input Wave Height (JONSWAP)'},'Location','southeast')
    F(ii)=getframe(gcf)
    pause(0.5)
    delete(p2)
    delete(p3)
    delete(p4)
    delete(p1)
    delete(p5) 
end
% create the video writer with 1 fps
  writerObj = VideoWriter('transect2_rd_int_15wr.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

