% ----------       refSpot v-0.7 ------------- % 
% first workable version of refSpot with the feature z-slice fitting
% calculate the intensity of individual spot (background substracted)
% refspot is designed to use microbeJ output, which help detect local
% maxima (i.e. fluorescent spot gaussian peak loacation, stored as a
% sub-pixel xy-coordinate pair), and then use Oufti output, which is a
% cell mesh, to generate binary mask (require the getBfData function from 
% Oufti) and the spot relative coordinates. Using the information of spot
% locations, the script will sum up 48 pixels (a 7x7 squre area) surround-

% ing that pixel and calculate average cell background intensity.
% Bkg-subtracted spot intensity will be stored in the Oufti cellList 
% stucture under each cell as individual rows in the maxtrix. The structure 
% field is called 'spot'.
% --- this version does not include the feature of finding maxZ --- %
% version: 0.7a
% -----------
% v0.6 - Adding universal file name feature
% v0.6b: exclude spots on the edges
% Including a new feature of counting MOI
% Required function: getBfData, projectToMesh
% Author: Zihao Yu
% Date created: 02/19/21 
% Last modified: 12/10/23

%% Main
clear
clc

load ouftiOut;
load microbeJ_max;
spotMax = test;
clear test;

load microbeJ_all;
spotAll = test;
clear test;

cd('signal/'); % could be signal or c2 which represent for the second/fluorescence channel
start = 1;
startA = 1;
freePhage = 0;
s2 = strel('disk',3);

for nframe = 1:spotMax(end).POSITION.position
    imSignal = imread(['red' num2str(nframe-1,'%04d') '.tif']); % whether to use nframe or nframe-1 depend on the actual beginning number, e.g., green00 or green01
    for ncell = 1:length(cellList.meshData{nframe})
        if cellList.meshData{nframe}{ncell}.mesh ~=0 
            cellStructure = cellList.meshData{nframe}{ncell};
            [cellList.meshData{nframe}{ncell},polymask,cellMask,cellCore] = getBfData(cellStructure); % modified from getextraData and previously used calculation.m code
            cellList.meshData{nframe}{ncell}.spot = []; % create the field to store spot information
            tempabs = [];
            temprel = [];
            allSpot = [];
            counter = 0;
            rawSpotInt = [];
            spotMask = zeros(512,512);
            for nspot = start:length(spotMax)
                position = spotMax(nspot).POSITION.position;
                if position == nframe
                    % -- get absolute coordinates of spots
                    spotX = round(spotMax(nspot).GAUSSIAN.x);
                    spotY = round(spotMax(nspot).GAUSSIAN.y);

                    if inROI(polymask,spotX,spotY)
                        if spotX > 3 && spotX < 510 && spotY >3 && spotY < 510
                            tempabs = [tempabs;spotX spotY];
                            rawSpotInt = [rawSpotInt sum(imSignal(spotY-3:spotY+3, spotX-3:spotX+3),'all')];
                            % -- get relative coordinates in Oufti mesh projection format
                            [l,d] = projectToMesh(spotX,spotY,cellStructure.mesh); % check more info by typing help projectToMesh
                            temprel = [temprel;l d];
                            counter = counter + 1;
                        end
                    else 
                        freePhage = freePhage + 1;
                    end
                elseif position ~= nframe
                    break
                end
            end
            nspotM = nspot;
            clear nspot;
            cd ('../red');
            c2Files = dir(pwd);
            filePrefix = extractBefore(c2Files(4).name,'1c'); % files(3) to skip system folders, '01c' to help set the position to get file prefix (apply to standard xy c t/z NIS TIFF export)
            channelNo = extractAfter(extractBefore(c2Files(4).name,'z'),'1c');
            for nspot = startA:length(spotAll)
                if spotAll(nspot).POSITION.slice > 5*(nframe-1) && spotAll(nspot).POSITION.slice < 5*nframe+1 % spots from all z slices but in the same frame
                    zSlice = double(spotAll(nspot).POSITION.slice - 5*(nframe-1));
                    imZ = imread([ filePrefix num2str(nframe,'%01d') 'c' channelNo 'z' num2str(zSlice,'%01d') '.tif']);
                    spotX = round(spotAll(nspot).GAUSSIAN.x);
                    spotY = round(spotAll(nspot).GAUSSIAN.y);                    
                    if inROI(polymask,spotX,spotY)  % spots from all z slices but in this cell only
                        if spotX > 3 && spotX < 510 && spotY >3 && spotY < 510
                            allSpot = [allSpot;spotX spotY zSlice sum(imZ(spotY-3:spotY+3, spotX-3:spotX+3),'all')];
                        end
                    end
                end
            end
            clear spotX;
            clear spotY;
            nspotA = nspot;
            clear nspot;
            spotZfit = {};
            for i = 1:counter
                spotZtemp = [];
                for j = 1:size(allSpot,1)
                    dist = sqrt(sum(allSpot(j,1:2)-tempabs(i,:)).^2);
                    if dist < 4
                        spotZtemp = [spotZtemp;allSpot(j,:)];
                    end
                end
                if ~isempty(spotZtemp)
                    [maxInt,spotIndex] = max(spotZtemp(:,4));
                    spotZfit{i} = spotZtemp(spotIndex,:);
                else
                    spotZfit{i} = [tempabs(i,:) 0 rawSpotInt(i)];
                end
                spotX = spotZfit{i}(1);
                spotY = spotZfit{i}(2);
                % -- round spots on the edges to create masks -- %
                if spotX > 507 
                    spotX = 507;
                elseif spotY > 507
                    spotY = 507;
                elseif spotX < 6
                    spotX = 6;
                elseif spotY < 6
                    spotY = 6;
                end
                spotMask(spotY-5:spotY+5, spotX-5:spotX+5) = 1;
                clear spotZtemp;
            end
            cellList.meshData{nframe}{ncell}.spot.location_maxProj = tempabs; % 02/16 update adding '-maxProj' to distinguish
            cellList.meshData{nframe}{ncell}.spot.meshProjection_maxProj = temprel;
            cellList.meshData{nframe}{ncell}.spot.number = counter;
            cellList.meshData{nframe}{ncell}.spot.location_all = spotZfit;
            spotInt = {};
            if ~isempty(cellList.meshData{nframe}{ncell}.spot)
                for actualSpot = 1:cellList.meshData{nframe}{ncell}.spot.number
                    cellBkgMask = cellMask - spotMask;
                    %                 figure; imshow(cellBkgMask);
                    imCellBkg = immultiply(cellBkgMask == 1, imSignal);
                    cellList.meshData{nframe}{ncell}.cellBkgInt = sum(imCellBkg(:));
                    cellBkgArea = sum(cellBkgMask(:));
                    spotBkgInt = (cellList.meshData{nframe}{ncell}.cellBkgInt/cellBkgArea)*49; % the number 49 represent for 7x7 squre area
                    spotInt{actualSpot} = spotZfit{actualSpot}(4) - 0;
                end
            end
            actualInt = cell2mat(spotInt);
            cellList.meshData{nframe}{ncell}.spot.intensity = actualInt;
        end
    end
    start = nspotM;
    startA = nspotA;
    cd ..;
end


close;
save(['cellInfo'],'cellList');

%% xProj localization
load cellInfo.mat;
xLoc = [];
count = 0;

for nFrame = 1:length(cellList.meshData);
    for nCell = 1:length(cellList.meshData{nFrame})
        if cellList.meshData{nFrame}{nCell}.mesh ~=0
            if cellList.meshData{nFrame}{nCell}.spot.number ==1
                spotLoc = cellList.meshData{nFrame}{nCell}.spot.meshProjection_maxProj;
                xLoc = [xLoc; spotLoc(1)./cellList.meshData{nFrame}{nCell}.length];
            elseif cellList.meshData{nFrame}{nCell}.spot.number >1
                if cellList.meshData{nFrame}{nCell}.spot.number  == 2
                    tempLoc = cellList.meshData{nFrame}{nCell}.spot.meshProjection_maxProj;
                    tempX = tempLoc./cellList.meshData{nFrame}{nCell}.length;
                    if (tempX(1,1)>0.1 && tempX(1,1)<0.8) || (tempX(2,1)>0.1 && tempX(2,1)<0.8)
                        count = count+1;
                    else
                        xLoc = [xLoc; tempX(:,1)];
                    end
                else
                    count = count+1;
                end
            end
        end
    end
end

cellRefX = abs(xLoc-0.5)*2;
polePP7 = (sum(cellRefX<0.1)+sum(cellRefX>0.8))/length(cellRefX);