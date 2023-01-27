clear all;

%-------------------
%Set up floodplain elevation points?
%------------------------

% w = dlmread('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/floodplain.xyz','|'); %detects the delimiter from the file and treats repeated white spaces as a single delimiter
% 
% disp(['Finished reading buffer point file at clock set: ', num2str(clock),'...']);
% 
% [proj] = geotiffinfo('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/wleopold3k_check_ARK_TMRA.tif');
% [fldpln_x, fldpln_y] = projfwd(proj, w(:,2), w(:,1));
% 
% %store as xyz
% fldplnbuffer(:,1) = fldpln_x;
% fldplnbuffer(:,2) = fldpln_y;
% fldplnbuffer(:,3) = w(:,3);
% 
% save('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/fldpln_buffer_pts.mat','fldplnbuffer');




%--------------------------
%Load Data
%---------------------------


%floodplain edge points?
fp = open('/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/downscaling/fldpln_buffer_pts.mat');

% load validation data. pull data and metadata from ascii raster
    %w3k: array
    
[w3k, ncols3k, nrows3k, xllcorner3k, yllcorner3k, cellsize3k, nodata] = ascii_reader(['/Users/gschuman/Documents/JPL_postdoc_work/LISFLOOD/RedArkansas/4Kostas/LISFLOOD_Truth/wleopold3k_check_ARK.asc']);

%load terrain data?
 
numfiles = 30; %time steps (output from LISFLOOD)

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

       
    %get indicies where sgwsl is REAL and w3k is NOT real
    %"exclude those 30 m cells that represent the bankfull channel in sub-grid scale according to Landsat"
     [f3k_r f3k_c] = find(sgwsl~=-9999 & w3k==-9999); % returns the row and column subscripts of each nonzero element 
     
    %%%% operate on floodplain %%%%
 
    
    %convert sgwsl indicies to spatial mid-cell
    %starting points?
    yul3k = yllcorner3k+(nrows3k*cellsize3k)+cellsize3k;
    xul3k = xllcorner3k;    
    
    f3k_xyz(:,1) = xul3k+((f3k_c*cellsize3k)-(cellsize3k)/2); %mid cell. x locations. convert to spatial
    f3k_xyz(:,2) = yul3k-((f3k_r*cellsize3k)+(cellsize3k)/2); %mid cell. y locations?


    %identify all the sgwsl indicies within range of the floodplain buffer
    %https://www.mathworks.com/help/stats/rangesearch.html
    %default value is 'exhaustive'.  software computes the distances from all X points to each Y point to find nearest neighbors.
    [id,Dis]= rangesearch( %Find all neighbors within specified distance
        fp.fldplnbuffer(:,1:2), %identify these points (floodplain edge?). fine resolution. 
        f3k_xyz(:,1:2), %within range of these points. coarse resolution
        cellsize3k*1.5, %range
        'Distance','Cityblock');
        
        
    disp(['Found all floodplain points for sgwsl ', num2str(i),' in high res file for each low res floodplain point at clock set: ', num2str(clock),'...']);

    fpbuff = fp.fldplnbuffer;

    %loop through each coarse sgwsl point and compute differences with neighbouring hi-res
    %"look for the closest 30 m cells whose elevations are lower than the simulated water surface elevation"
    for k = 1:length(id) 
    
        %get this floodplain (fine) point info
        getxyzid = fpbuff(
            id{k,1}, %id of best matching fine point
            1:3); %get xyz's for this set
            
        getxyzid(:,4) = id{k,1}; %and the id
        
        %get the low-res WSL (z) value 
        z_sgwsl = sgwsl(f3k_r(k),f3k_c(k)); 
        
        %loop on each matching hi-res point
        for kk =  1:length(getxyzid);
            
            %hi-res WSL less than low-res
            if getxyzid(kk,3)<z_sgwsl
            
                %change this hi-res point 
                fpbuff(
                    getxyzid(kk,4), %hi-res ID
                    4 %new column
                    ) = z_sgwsl-getxyzid(kk,3); %low-res WSE - hi-res WSE
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


