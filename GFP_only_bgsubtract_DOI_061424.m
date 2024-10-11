%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================

% -- FIRST, combine all FOVs horizonally in ImageJ/Fiji -- %
% -- Update the number of conditions and FOVs in the script -- %
function GFP_only_bgsubtract_DOI_061424(filename_base)

%clear all
%close all

% load the image tif files data
filepath='/Users/tinakim/Dropbox/MATLAB/SOPHIA/Data/';
filefolder='10142022_TREGFP_quantDOI_reanalysis_061424/';
filename=strcat(filename_base,'.tif');

imsegm=imread(strcat(filepath,filefolder,filename));
imsegm=double(imsegm);


% This assumes FOVs are combined in 1 row x nFOV columns.
%nFOV=20; % e.g., 10 -DOI and 10 +DOI FOVs in 1 row = 20
nFOV=8;
nCON=1; % e.g., only +light in one row
imx=size(imsegm,2)/nFOV; % width
imy=size(imsegm,1)/nCON; % height
FOV_means=zeros(nCON,nFOV);
FOV_cellcount=zeros(nCON,nFOV);


% Segmentation by adaptive thresholding
prm.method = 'adth';
prm.adth.th = 0.05;

[cellbw1,wat,imsegmout,prmout] = cellsegm.segmct(imsegm,0.02,0.3,'prm',prm);
cellsegm.show(imsegmout,1);title('Raw image');axis off;
cellsegm.show(cellbw1,1);title('Cell segmentation by ADTH');axis off;

% improving the results by splitting of cells
splitth = 1;
plane = 1;
% cells above this threshold are split (all cells here)
n = 100;
h = [0.5 0.5 1.5];
cellbw2 = cellsegm.splitcells(cellbw1,splitth,n,h);
cellsegm.show(cellbw2,2);title('Cell segmentation by ADTH with splitting');axis off;
%cellbw2=cellbw1;

% create cellmasks
[L,num_cells]=bwlabel(cellbw2);
segcentroids=zeros(num_cells,2);
meanfluo=zeros(num_cells,1);

for a=1:num_cells
    [tempx,tempy]=find(L==a);
    fluo=[]; % calculate each cell's fluorescence (e.g. scFLARE2 staining)
    cellmasktemp=zeros(size(cellbw2,1),size(cellbw2,2));
    for b=1:length(tempx)
        cellmasktemp(tempx(b),tempy(b))=1;
        fluo(b)=imsegm(tempx(b),tempy(b));
    end
    temp=regionprops(cellmasktemp);
    segcentroids(a,:)=temp.Centroid;
    meanfluo(a)=mean(fluo);
end

No_DOI_inds=find(segcentroids(:,1)<=imx*nFOV/2);
DOI_inds=find(segcentroids(:,1)>imx*nFOV/2);

% save data 
savename=strcat(filename_base,'_bgsub_fluo.mat');
save(strcat(filepath,filefolder,'results/',savename),'num_cells','meanfluo','No_DOI_inds','DOI_inds');

%% Calculate the average ratios per FOV

for z=1:nCON % number of rows in combined image
    temp=find(segcentroids(:,2)>=(z-1)*imy+1 & segcentroids(:,2)<(z)*imy); %finds cells corresponding to current row
    tempx=segcentroids(temp,1); %gets x position of each cell in current row
    
    imsegm_row=imsegm((z-1)*imy+1:(z)*imy,:);
    
    NoDIO=[];
    DOI=[];
    
    for a=1:nFOV % number of columns in combined image
        
        tempa=find(tempx>=(a-1)*imx+1 & tempx<(a)*imx); %finds cells in current row for current column
        currmask=cellbw2((z-1)*imy+1:z*imy,(a-1)*imx+1:a*imx);
        currim=imsegm_row(:,(a-1)*imx+1:a*imx);
        currim_bg=mean(currim(~currmask));
        %currim_bg=prctile(currim(~currmask),5);
        %currim_bg=0;
        
        cellfluo=meanfluo(temp(tempa));
        
        if a<(nFOV/2)+1
            NoDIO=[NoDIO; cellfluo];
        else
            DOI=[DOI; cellfluo];
        end
        
        FOV_means(z,a)=mean(cellfluo);

        FOV_cellcount(z,a)=length(tempa);
        
    end
end

savename=strcat(filename_base,'_FOV_bgsub_fluo.mat');
save(strcat(filepath,filefolder,'results/',savename),'FOV_means','FOV_cellcount','DOI','NoDIO');

display('done!')

end
