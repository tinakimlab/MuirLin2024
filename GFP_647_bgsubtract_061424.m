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
function GFP_647_bgsubtract_061424(filename_base)

%clear all
%close all

% load the image tif files data
filepath='/Users/tinakim/Dropbox/MATLAB/SOPHIA/Data/';
filefolder='10142022_all_quantDOI_reanalysis_061424/';
filename_static=strcat(filename_base,'_647.tif');

imsegm=imread(strcat(filepath,filefolder,filename_static));
imsegm=double(imsegm);
filename_TRE=strcat(filename_base,'.tif');
imsegm_TRE=imread(strcat(filepath,filefolder,filename_TRE));
imsegm_TRE=double(imsegm_TRE);


% This assumes FOVs are combined in 1 row x nFOV columns.
%nFOV=20; % e.g., 10 -calcium and 10 +calcium FOVs in 1 row = 20
nFOV=8;
nCON=1; % e.g., only +Biotin in 1 row
imx=size(imsegm,2)/nFOV; % width
imy=size(imsegm,1)/nCON; % height
FOV_ratios=zeros(nCON,nFOV);


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
meanfluo_TRE=zeros(num_cells,1);
meanratio=zeros(num_cells,1);
for a=1:num_cells
    [tempx,tempy]=find(L==a);
    fluo=[]; % calculate each cell's tTA fluorescence
    fluo_TRE=[]; % calculate each cell's TRE fluroescence
    cellmasktemp=zeros(size(cellbw2,1),size(cellbw2,2));
    for b=1:length(tempx)
        cellmasktemp(tempx(b),tempy(b))=1;
        fluo(b)=imsegm(tempx(b),tempy(b));
        fluo_TRE(b)=imsegm_TRE(tempx(b),tempy(b));
    end
    temp=regionprops(cellmasktemp);
    segcentroids(a,:)=temp.Centroid;
    meanfluo(a)=mean(fluo);
    meanfluo_TRE(a)=mean(fluo_TRE);
    meanratio(a)=mean(fluo_TRE)./mean(fluo); % calculate average TRE/tTA ratio
end

no_ca_inds=find(segcentroids(:,1)<=imx*nFOV/2);
ca_inds=find(segcentroids(:,1)>imx*nFOV/2);

% save data 
savename=strcat(filename_base,'_bgsub_ratios_test.mat');
save(strcat(filepath,filefolder,'results/',savename),'num_cells','meanfluo','meanfluo_TRE','meanratio','no_ca_inds','ca_inds');

%% Calculate the average ratios per FOV

for z=1:nCON % number of rows in combined image
    temp=find(segcentroids(:,2)>=(z-1)*imy+1 & segcentroids(:,2)<(z)*imy); %finds cells corresponding to current row
    tempx=segcentroids(temp,1); %gets x position of each cell in current row
    
    imsegm_row=imsegm_TRE((z-1)*imy+1:(z)*imy,:);
    
    Light_NoCa=[];
    Light_Ca=[];
    
    for a=1:nFOV % number of columns in combined image
        
        tempa=find(tempx>=(a-1)*imx+1 & tempx<(a)*imx); %finds cells in current row for current column
        currmask=cellbw2((z-1)*imy+1:z*imy,(a-1)*imx+1:a*imx);
        currim=imsegm_row(:,(a-1)*imx+1:a*imx);
        currim_bg=mean(currim(~currmask));
        %currim_bg=prctile(currim(~currmask),5);
        %currim_bg=0;
        
        tTA=meanfluo(temp(tempa));
        TRE=meanfluo_TRE(temp(tempa))-currim_bg;
        
        if a<(nFOV/2)+1
            Light_NoCa=[Light_NoCa; TRE./tTA];
        else
            Light_Ca=[Light_Ca; TRE./tTA];
        end
        
        FOV_ratios(z,a)=mean(TRE)/mean(tTA);

        % calculate mean cellmasks
        %FOV_ratios(z,a)=mean(currim(currmask))-mean(currim_bg);
        %FOV_cellcount(z,a)=length(tempa);

%         figure(100);
%         imagesc(currim);
%         disp(mean(currim(currmask)));
%         disp(mean(currim_bg));
%         pause;
        
    end
end

savename=strcat(filename_base,'_FOV_bgsub_ratios_test.mat');
save(strcat(filepath,filefolder,'results/',savename),'FOV_ratios','Light_Ca','Light_NoCa');

display('done!')

end
