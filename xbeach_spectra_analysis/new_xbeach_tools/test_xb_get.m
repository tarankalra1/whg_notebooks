%Code to generate the wave heights
%% Spectral analysis of waterlevel inputs
close all
clear all

%ncfile='xboutput.nc';

ncfile='xboutputE.nc';

ncinfo=ncinfo(ncfile);
t=ncread(ncfile,'globaltime');
x=ncread(ncfile,'globalx');
zs=ncread(ncfile,'zs');
zs=squeeze(zs)';

xbo = xb_get_spectrum(zs)

hh = xbo.data ;
Hm0=hh(7).value ;
Snn=hh(4).value ;
save('xb_get_run_Emily.mat','x','Snn','Hm0');

