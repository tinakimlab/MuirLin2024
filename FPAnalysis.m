%% FP Script
%Describes analysis performed in Muir&Lin et al 2024 (figure 1B-C)

%Import Files
filepath=['C:\Users\Jessie\Desktop\PFC-DOI-scFlare\PsychLight\DelayedInjection\']
files=dir(filepath)
filename={}
m=1
for f=3:length(files)
    filename{m}=files(f).name
    m=m+1
end

salInj=%Insert timestamp of Injection for each recording
DOIInj=%Insert timestamp of Injection for each recording
raw=[];
con=[];
dFF=[];
time=[]
%Calculate dff and mean fluorescence in pre-determined time bins (Figure
%1C) for each recording and save in array
for n=1:length(filename)
    data=[filepath,char(filename(n))]
    M=table2array(readtable(data));
    raw=M(2:2:144006,3);
    con=M(1:2:144006,3);
    time=(M(1:2:144006,1)*10^-3)/60;
    dFF(n,:)=deltaFF(raw,con); 
%Calculate mean fluorescence in time bins- adjust time bins for specific
%experiments
    [a,salInd]=min(abs(time-salInj))
    [a,DOIInd]=min(abs(time-DOIInj))
    base(n)=mean(dFF(1:salInd(n)),1)
    sal(n)=mean(dFF(salInd(n):DOIInd(n)),1)
    doi(n,:)=[mean(dFF(DOIInd(n):DOIInd(n)+18002),1);mean(dFF(DOIInd(n)+18003:DOIInd(n)+36004),1);
        mean(dFF(DOIInd(n)+36005:DOIInd(n)+54006),1);mean(dFF(DOIInd(n)+54006:DOIInd(n)+72008),1);
        mean(dFF(DOIInd(n)+72008:DOIInd(n)+90010),1);mean(dFF(DOIInd(n)+90011:DOIInd(n)+108012),1);
        mean(dFF(DOIInd(n)+108013:DOIInd(n)+126014),1);mean(dFF(DOIInd(n)+126014:end),1)]

 end


% Plot average traces with SEM shaded (Figure 1B, S1B,C)
avedFFT=zscore(mean(dFF,1));
sem=std(avedFFT,0,2)/sqrt(5)
dff1 =smooth((avedFFT-sem),20);
dff2 = smooth((avedFFT+sem),20);
inBetween = [dff1; flipud(dff2)];
x(:,1)=[1:length(time)];
x2=[x(:,1); flipud(x(:,1))];
figure
fill(x2, inBetween, 'k', 'LineStyle', 'none','FaceAlpha',0.5);

hold on
plot(time/60, smooth(avedFFT,20), 'k', 'LineWidth', 0.5);
plot(x, dff1,'b',x,dff2, 'b', 'LineWidth', 0.5);
ylim([-3 3])
hold off
 print('YourEPSFile','-desp','-vector')

