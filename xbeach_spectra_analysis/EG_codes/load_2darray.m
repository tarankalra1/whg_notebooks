%clear all ; close all ; clc ;

load('100yr_degraded_transect4_Hm.mat')
% the time of chunks for Hm0 
for t = 1:length(Hm0)
%    data_hf(i) = Hm0_hf{i};
%    data_lf(i) = Hm0_lf{i};
  data0       = Hm0{t};
  data1(t,:)  = data0 ;
%    time1(i)   = window*(i); 
    %time1(1,i}     = window*t_int{1,1};  
end  
%my_array = reshape(cell2mat(Hm0), [312, 1245]);
plot(data1(:,203))

save('100yr_degraded_transect4_processed.mat','data1')