clear;
file=['C:\Users\Jessie\F2dbdata.csv']
M=table2array(readtable(file));

l=607 %Frames to analyze (10 min)-change based on frequency
r=891 %length of recording
Mnorm=[];
SDratio=[]
DDratio=[]
filter=0.007%0.008
%Detrend data, calculate difference between saline and DOI (SDratio) and
%30min DOI and 0min DOI (DDratio)
for t=1:length(M(:,1))
    Mnorm(t,:)=highpass(M(t,:),filter,1);
    SDratio(t,:)=Mnorm(t,r+1:r+l)-Mnorm(t,1:l);
    DDratio(t,:)=Mnorm(t,r*2+1:r*2+l)-Mnorm(t,1:l);
end

%Calculate Average difference for both periods and sort
ave=mean(SDratio,2);
aveDD=mean(DDratio,2)
[ave,sortIdx] = sort(ave,'descend');
aveDD=aveDD(sortIdx)

doiCell=[]
salineCell=[]
neutralCell=[]
thresh=max(ave(2:end))*0.05;%Threshold
%Use thresh to sort cells in 0-10 period
for c=1:length(ave)
    if ave(c)>=thresh
        doiCell=[doiCell,c];
    elseif ave(c)<=-thresh
        salineCell=[salineCell,c];
    else 
        neutralCell=[neutralCell,c];
    end
end

doi30Cell=[]
saline30Cell=[]
neutral30Cell=[]
%Use same thresh to sort cells in 30-40 period
for c30=1:length(aveDD)
    if aveDD(c30)>=thresh
        doi30Cell=[doi30Cell,c30];
    elseif aveDD(c30)<=-thresh
        saline30Cell=[saline30Cell,c30];
    else 
        neutral30Cell=[neutral30Cell,c30];
    end
end
%Find intersection
doubleAct=intersect(doi30Cell,doiCell);
