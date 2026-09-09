%% this is a parent function that will read in a given path to a stacked dicom (nx, ny, nz*b) and output a set of maps following a certain algorithm
% input: path to dicom
% assumes bvalues, averages, and the b-value cutoff (where you want to split between perfusion and diffusion)
% also attempts automatic segmentation


%b=[0,50,100,150,200,250,300,500,700,900,1000]
% guessing: 4 averages for each nonzero bvalue, no averages for b=0? 
% this means (10*4 + 1) * n_slices so for 1025, it's 25 slices... 


% Mira Liu Sept 8 2026

function IVIM_maps = RunIVIM_Algos(dicompath,dicom)

    DicomStack = dicomread([dicompath,dicom]);
    DicomStack = squeeze(DicomStack);
    [nx,ny,nz] = size(DicomStack);

    % hardcodign now, but can make it an input if you want to change it. 
    n_averages = 4;
    bvalues = [0,50,100,150,200,250,300,500,700,900,1000];
    Images_per_slice = n_averages.*(length(bvalues)-1)+ 1;
    bval_cutoff_idx = 6; % so 250 is the first data point considered perfusion  

    n_slices = nz/Images_per_slice;

   
    % sort the data into a nx by ny by slice by b-value matrix
    temp_grouped = reshape(DicomStack, [nx,ny, Images_per_slice, n_slices]);
    slices_grouped  = permute(temp_grouped, [1,2,4,3]);
    b0 = slices_grouped(:,:,:,1);
    b_nonzero = slices_grouped(:,:,:,2:end);
    b_nonzero_split = reshape(b_nonzero, [nx, ny, n_slices, n_averages, length(bvalues)-1]);
    b_nonzero_averaged = mean(b_nonzero_split,4);
    b_nonzero_averaged = squeeze(b_nonzero_averaged);
    DicomStack = cat(4, b0, b_nonzero_averaged);


    % now run the analysis to create a set of maps (f, D, Dstar) for each slice.
    IVIM_maps = zeros(nx,ny,n_slices,6);
    for slice = 1:n_slices
        for i = 1:nx
            for j = 1:ny
                if DicomStack(i,j,slice,1) > 150 % attempt at brain mask-ish
                    signal = double(squeeze(DicomStack(i,j,slice,:)));
                    Output = Algorithm6(bvalues, signal,bval_cutoff_idx); %here you can change algorithms if you want, Algo6 is the one we decided on
                    IVIM_maps(i,j,slice,1) = Output.f;
                    IVIM_maps(i,j,slice,2) = Output.D;
                    IVIM_maps(i,j,slice,3) = Output.Dstar;
                    IVIM_maps(i,j,slice,4) = Output.SSE;
                    IVIM_maps(i,j,slice,5) = Output.rsq;
                    IVIM_maps(i,j,slice,6) = Output.adj_rsq;
                end
            end
        end
    end


    % now also make a brain mask from b=1000
    b_ind = 11; % the bvalue index for b=1000
    brainMask = segmentBrain_IVIM(double(DicomStack),b_ind);

    % now save them
    save ([dicompath, 'brainmask'], 'brainMask')
    save([dicompath,'IVIM_Map'], 'IVIM_maps')


end




function [brainMask] = segmentBrain_IVIM_map(imgstack,b_ind)
%% Segment IVIM images
% Adopted from Yong Jeong SegmentBrain_V2
% ML 07/14/2022

[nr,nc,ns,nb,ndgo] = size(imgstack);
imgstack = squeeze(imgstack(:,:,:,b_ind,:));
meanMat = mean(imgstack,4);
sigmaMat = std(imgstack,0,4);

cvMat = (sigmaMat./meanMat) < .13; % coefficient of variation
brainMask = zeros(nr,nc,ns);

% Get biggest contiguous 3d mask
brainMask = getLCMask(cvMat);
% Find the biggest area slc and centroid
biggestArea = 0;
biggestCentLoc = [];
for ii = 1:ns
    tmpmask = getLCMask(brainMask(:,:,ii));
    s = regionprops(tmpmask,'Centroid');
    currarea = regionprops(tmpmask,'area');
    if ~isempty(currarea)
        if currarea.Area > biggestArea
            biggestArea = currarea.Area;
            biggestCentLoc = round(s.Centroid);
            biggestSlc = ii;
            se = strel('disk',1,4);
            tmpmask = imdilate(tmpmask,se);
            tmpmask = imfill(tmpmask,'holes');
            biggestSlcMask = tmpmask;
        end
    end
end
 
% Remove regions that do not overlap with biggest slc mask
for ii = 1:ns
    tmpmask = brainMask(:,:,ii);
    s = bwconncomp(tmpmask);
    if iscell(s.PixelIdxList)
        numregion = length(s.PixelIdxList);
        regionidxlist = s.PixelIdxList;
    else
        numregion = 1;
        regionidxlist = {s.PixelIdxList};
    end
    for jj = 1:numregion
        if isempty(intersect(regionidxlist{jj},find(biggestSlcMask)))
            tmpmask(regionidxlist{jj}) = 0;
        end
        brainMask(:,:,ii) = tmpmask;
    end
end        
    
% biggestSlcMask = repmat(biggestSlcMask,[1 1 ns]);
 brainMask = biggestSlcMask.*brainMask;

 

% Erode mask to remove appendages...
se = strel('disk',1,4);
brainMask = imerode(brainMask,se);
brainMask = imerode(brainMask,se);

% Get biggest contiguous 2d mask for each slice
for ii = 1:ns
    tmpmask = getLCMask(brainMask(:,:,ii));
    brainMask(:,:,ii) = tmpmask;
end

% Dilate the mask and fill in holes

se = strel('disk',1,4);
brainMask = imdilate(brainMask,se);
for ii = 1:ns
    brainMask(:,:,ii) = imfill(brainMask(:,:,ii),'holes');
end
%}

% Dilate the mask
%se = strel('disk',1,4);
%brainMask = imdilate(brainMask,se);

% Remove masks that do not overlap with biggest area centroid (most
% likely stuff that's off the brain
for ii = 1:ns
    if brainMask(biggestCentLoc(2),biggestCentLoc(1),ii) == 0
        brainMask(:,:,ii) = zeros(nr,nc);
    end
end

% Get biggest contiguous 3d mask
brainMask = getLCMask(brainMask);

 %}
end

