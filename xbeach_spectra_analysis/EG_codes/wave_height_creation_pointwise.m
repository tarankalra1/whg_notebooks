%Code to generate the wave heights
%% Spectral analysis of waterlevel inputs
close all
clear all
clc ;
load('data_WSE.mat')
t=time;
zs = wse_point1 ; 

%wave record 15 minute time period is 900sec.
window_size = 15; 
window = window_size*60; 

delta_t=floor(length(zs)/window) ; 

ts=1 ;
for i=1:delta_t
   % n1(i)=ts;
    n=ts;
    ts=ts+window-1; 
    zs_i=zs(n:ts);
    zs_int{i}=zs_i;
     
    t_i=t(n:ts);
    t_int{i}=t_i;
end


index_j=length(zs_int);
k=1 ; 
for j=1:index_j
    zs_i2=zs_int{j};
    
    zs_i2_p = zs_i2'; 
    t_i2=t_int{j};
    
    [f,Snn,A,ff,Af] = comp_spec( detrend(zs_i2_p(:,1)), 1/(t_i2(2)-t_i2(1))  );
    df = f(2)-f(1);

    % calculate variance below (IG) and above (short waves) fsplit
    fsplit = 0.05;
    
    m0hf(k) = sum(Snn(f >= fsplit)*df);% short wave energy density variance
    m0lf(k) = sum(Snn(f <  fsplit)*df);%  long wave energy density variance
                
    Hm0_hfi = 4*sqrt(m0hf);
    Hm0_lfi = 4*sqrt(m0lf);
    Hm0_i = sqrt(Hm0_hfi.^2+Hm0_lfi.^2);
       
    Hm0_hf{j}=Hm0_hfi;
    Hm0_lf{j}=Hm0_lfi;  
    Hm0{j}=Hm0_i;
end 

for i = 1:length(Hm0);
    data_hf(i) = Hm0_hf{i};
    data_lf(i) = Hm0_lf{i};
    data0(i)   = Hm0{i};
    time1(i)   = window*(i); 
    %time1(1,i}     = window*t_int{1,1}; 
end 
figure(1) 
plot(time1/(3600), data0,'k')
hold on; 
plot(time1/(3600), data_lf, 'r')
hold on ; 
plot(time1/(3600), data_hf, 'g')
xlim([0, 80])
legend('combined','lf','hf')
%title(['point 1 at num2str(window} min'])
%title(sprintf('point 1 at %d', window_size, 'min'));
title(sprintf('Point1 at %d minutes', window_size));
%title(sprintf('Point 1 at %d minutes', window_size));
saveas(gcf, sprintf('Point1_%d_minutes.png', window_size), 'png');

%print(gcf, '-dpng', f"data_{window}.png")