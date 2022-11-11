clear all;
% w = dlmread('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/floodplain.xyz','|');
% 
% disp(['Finished reading buffer point file at clock set: ', num2str(clock),'...']);
% 
% [proj] = geotiffinfo('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/wleopold3k_check_ARK_TMRA.tif');
% [fldpln_x, fldpln_y] = projfwd(proj, w(:,2), w(:,1));
% 
% fldplnbuffer(:,1) = fldpln_x;
% fldplnbuffer(:,2) = fldpln_y;
% fldplnbuffer(:,3) = w(:,3);
% 
% save('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/fldpln_buffer_pts.mat','fldplnbuffer');

fp = open('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/fldpln_buffer_pts.mat');
[w3k, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/4Kostas/LISFLOOD_Truth/wleopold3k_check_ARK.asc']);

    
numfiles = 30;

for j = 1:numfiles;
i = j+359;
    if i < 10;
[sgwsl, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/elevs/ark-000', num2str(i),'.elev']);
% [sgd] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/wd/ark-000', num2str(i),'.wd']);
    elseif i >= 10 && i <100;
[sgwsl, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/elevs/ark-00', num2str(i),'.elev']);
% [sgd] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/wd/ark-00', num2str(i),'.wd']);
    elseif i >= 100 && i <1000;
[sgwsl, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/elevs/ark-0', num2str(i),'.elev']);
% [sgd] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/wd/ark-0', num2str(i),'.wd']);
    else
[sgwsl, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/elevs/ark-', num2str(i),'.elev']);
% [sgd] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/wd/ark-', num2str(i),'.wd']);
    end

    
 [f3k_r f3k_c] = find(sgwsl~=-9999 & w3k==-9999);
 
    %%%% operate on floodplain %%%%
yul3k = yllcorner3k+(nrows3k*cellsize3k)+cellsize3k;
xul3k = xllcorner3k;

f3k_xyz(:,1) = xul3k+((f3k_c*cellsize3k)-(cellsize3k)/2); %mid cell
f3k_xyz(:,2) = yul3k-((f3k_r*cellsize3k)+(cellsize3k)/2); %mid cell

[id,Dis]= rangesearch(fp.fldplnbuffer(:,1:2),f3k_xyz(:,1:2),cellsize3k*1.5,'Distance','Cityblock');
disp(['Found all floodplain points for sgwsl ', num2str(i),' in high res file for each low res floodplain point at clock set: ', num2str(clock),'...']);

fpbuff = fp.fldplnbuffer;

for k = 1:length(id)
getxyzid = fpbuff(id{k,1},1:3);
getxyzid(:,4) = id{k,1};
z_sgwsl = sgwsl(f3k_r(k),f3k_c(k));
for kk =  1:length(getxyzid);
    if getxyzid(kk,3)<z_sgwsl
        fpbuff(getxyzid(kk,4),4) = z_sgwsl-getxyzid(kk,3);
    else
    end
end
clear z_sgwsl;
end 

disp(['Downscaled floodplain for simulation timestep ', num2str(i),' at clock set: ', num2str(clock),'...']);

filename = ['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/fullyearcycle/fldpln_sim_', num2str(i),'.mat'];
save(filename, 'fpbuff', '-ascii','-double');

disp(['Written file at clock set: ', num2str(clock),'...']);
end


