clear all;
ch_pts = open('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/ch_buffer_at_fldpln_pts.mat');
fp = open('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/fldpln_buffer_pts.mat');
[w3k, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/4Kostas/LISFLOOD_Truth/wleopold3k_check_ARK.asc']);

disp(['Finished reading buffer point file at clock set: ', num2str(clock),'...']);

%%%get coordinates of channel

w = ch_pts.w2;
xchan = w(:,1);
ychan = w(:,2);
indchan(1:length(w(:,1)),1) = (1:1:length(w(:,1)));

disp(['Finished reading coordinates of buffer point file at clock set: ', num2str(clock),'...']);

    %%%% operate on channel %%%%
    

%%% extract water level from xs mid point and put same water level on each xs point in xsallpts_chan

%%% get closest point per cross-section only and spread on other points with
%%% same xs cat number

numfiles = 370;

t = zeros(1,numfiles);

for j = 341:numfiles;
    tic
i = j+359;
    if i < 10;
[sgwsl, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/output.truth/ark-000', num2str(i),'.elev']);
[sgd] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/output.truth/ark-000', num2str(i),'.wd']);
    elseif i >= 10 && i <100;
[sgwsl, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/output.truth/ark-00', num2str(i),'.elev']);
[sgd] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/output.truth/ark-00', num2str(i),'.wd']);
    elseif i >= 100 && i <1000;
[sgwsl, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/output.truth/ark-0', num2str(i),'.elev']);
[sgd] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/output.truth/ark-0', num2str(i),'.wd']);
    else
[sgwsl, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/output.truth/ark-', num2str(i),'.elev']);
[sgd] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/output.truth/ark-', num2str(i),'.wd']);
    end

disp(['Starting downscaling of model simulation ', num2str(j),'/', num2str(numfiles),'...']);
    
xs = csvread('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/transects_ARK.csv',1,0);

disp(['Finished reading xs point file and model files at clock set: ', num2str(clock),'...']);
    
%%%% interpolate wsl for tinwsl coordinates in channel

locsID3k(:,1) = ceil((xs(:,4) - xllcorner3k)/cellsize3k);
locsID3k(:,2) = nrows3k - floor((xs(:,5) - yllcorner3k)/cellsize3k);

[C,ia,ic] = unique(xs(:,2),'rows');
%%% 
for nv = 1:length(ia(:,1));
f_id = find(xs(:,2) == xs(ia(nv,1),2));
for nvv = 1:length(f_id);
elev3k(nvv,1) = sgwsl(locsID3k(f_id(nvv,1),2),locsID3k(f_id(nvv,1),1)); %row = y, column = x
bathy3k(nvv,1) = sgwsl(locsID3k(f_id(nvv,1),2),locsID3k(f_id(nvv,1),1))-sgd(locsID3k(f_id(nvv,1),2),locsID3k(f_id(nvv,1),1)); %row = y, column = x
end
fmax = find(elev3k == max(elev3k));
xs(f_id,6) = max(elev3k);
xs(f_id,7) = bathy3k(fmax(1,1));
clear f_id elev3k bathy3k fmax;
end
clear nv;

disp(['Finished reading off 3km modeled values at clock set: ', num2str(clock),'...']);

%%%% interpolate wsl for tinwsl coordinates in channel

disp(['Interpolating bathy and WSL now... at clock set: ', num2str(clock),'...']);

xs_4interp = xs;
fnan = find(xs_4interp(:,6)==-9999);
xs_4interp(fnan,:) = []; 
clear fnan;

tin_chan(:,1) = griddata(xs_4interp(:,4),xs_4interp(:,5),xs_4interp(:,7),xchan,ychan);
tin_chan(:,2) = griddata(xs_4interp(:,4),xs_4interp(:,5),xs_4interp(:,6),xchan,ychan);

w(:,3) = tin_chan(:,1);
w(:,4) = tin_chan(:,2)-tin_chan(:,1);
fnan = isnan(w(:,3));
f4nan = find(fnan==1);
w(f4nan,3) = 0;
w(f4nan,4) = 0;
clear xs_4interp f4nan fnan 

disp(['Finished operation on channel at clock set: ', num2str(clock),'...']);

    %%%% Finished operation on channel %%%%

[f3k_r f3k_c] = find(sgwsl~=-9999 & w3k==-9999);
     
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
  %%%% Finished operation on floodplain %%%%
  
disp(['Downscaled channel and floodplain for simulation timestep ', num2str(i),' at clock set: ', num2str(clock),'...']);

%%%% Merge channel and floodplain %%%%%       

[C_x,indch_x,indfldpln_x] = intersect(w(:,1),fpbuff(:,1)); %% coordinates are so precise that there is only one exact match in x or y

fpbuff(indfldpln_x,:) = w(indch_x,1:4);

disp(['Merged channel and floodplain for simulation timestep ', num2str(i),' at clock set: ', num2str(clock),'...']);

filename = ['/Volumes/LaCie/Downscaling_ARK/fullyearcycle/tstep_', num2str(i),'.mat'];

save(filename,'fpbuff');

clear C_x Dis f3k_c f3k_r f3k_xyz fpbuff getxyzid ia ic id iii indch_x indfldpln_x k kk locsID3k nvv sgd sgwsl tin_chan xs filename;   

t(j) = toc;

disp(['Finished downscaling of model simulation ', num2str(j),'/', num2str(numfiles),'...']);

end
