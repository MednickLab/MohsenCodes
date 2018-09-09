% Developed by Mohsen Naji on Aug 4, 2017 for Sleep and Cognition lab, UCR
% Modified on Apr 18, 2018
% 1) browse edf; 2) browse stage file 3) confirm the freq ranges
% I assume the data is re-referenced
% _v2 for edf files including equally-sampled data 
% Ignore the output files for E1/E2, CHIN/EMG and ECG channels

function myeegpwr_v3()
ctu=0;
disp('This program assumes edf files contain re-referenced data(to mastoids or earlobes)');
while ctu==0
    [FileName,PathName] = uigetfile([pwd '/*.edf'],'Select the sleep edf file');
    [FileNameM,PathNameM] = uigetfile([pwd '/*.mat'],'Select the marker mat file');
    load([PathNameM FileNameM]);
    disp([FileNameM ' is loaded as sleep marker']);
    disp([FileName ' is loaded as the edf file ']);
    prompt='Did you select correctly? Y/N ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'Y';
    end
    if (str=='y') || (str=='Y')
        ctu=1;
    end
end
%%%
prompt = {'Delta','Theta','Alpha','Slow Sigma','Fast Sigma'};
dlg_title = 'Frequency definition';
num_lines = 1;
defaultans = {'[0.5 4]','[4 8]','[8 13]','[9 11]','[12 15]'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
f_delta=str2num(answer{1,1});
f_theta=str2num(answer{2,1});
f_alpha=str2num(answer{3,1});
f_slowsigma=str2num(answer{4,1});
f_fastsigma=str2num(answer{5,1});
%%%

[X,Channels,fs]=edf2mat_samefs([PathName FileName]);

wins=stageData.win;
mrk=stageData.stages;
mrk(end)=[];
t=find((mrk(2:end)-mrk(1:end-1))~=0);
smp=[0;t*30*fs;size(mrk,1)*30*fs]';
bnd=zeros(length(smp)-1,2);
for i=1:length(smp)-1
    bnd(i,:)=[smp(i)+1 smp(i+1)]; % beginning and end of each bout
end
bnd(end,2)=size(X,2);
Stage=mrk([1;t+1]); % bout sleep stage
duration=(bnd(:,2)-bnd(:,1)+1)/(fs*60); % bout duration in minute
P=cell(length(duration)-1,length(Channels));
spins={'\\','|','/','-'};
reverseStr = 'Power spectrum estimation ';
e1ch=zeros(length(Channels),1);
for i=1:length(Channels)
    e1ch(i)=(Channels{1,i}(1)=='E' && Channels{1,i}(2)=='1');
end
e2ch=zeros(length(Channels),1);
for i=1:length(Channels)
    e2ch(i)=(Channels{1,i}(1)=='E' && Channels{1,i}(2)=='2');
end
emgch=zeros(length(Channels),1);
for i=1:length(Channels)
    emgch(i)=(Channels{1,i}(3)=='G' && Channels{1,i}(2)=='M'&& Channels{1,i}(1)=='E');
end
rj=[find(e1ch==1) find(e2ch==1) find(emgch==1)];
X(rj,:)=[];
Channels(rj)=[];

[b2,a2]=butter(4,0.16/(fs/2),'high');
[b3,a3]=butter(4,35/(fs/2),'low');
X_raw=filtfilt(b3,a3,X');
X_raw=filtfilt(b2,a2,X_raw);
X_raw=X_raw';
for i=1:length(duration)
    si=mod(i,4); 
    if si==0
        si=4;
    end
    msg = (spins{si});
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, 1);
    
    x=X_raw(:,bnd(i,1):bnd(i,2));
%     rj=[];
%     for ii=1:size(X_raw,1)
%         rj=[rj find(abs(x(ii,:))>320)];
%     end
%     irj=[];
%     for rr=1:length(rj)
%         irj=[irj max(1,rj(rr)-2*fs):min(size(x,2),rj(rr)+2*fs)];
%     end
%     x(:,irj)=[];
    p2=2.^[8:13];[~,kk]=min(abs(4*fs-p2));hw=p2(kk);    
    [pxx,f]=pwelch(detrend(x'),hanning(hw),hw/2,0:fs/hw:(fs/2),fs);
    for j=1:length(Channels)
        P{i,j}=pxx(:,j);
    end
end
Stage(end)=[];
duration(end)=[];
fprintf('\n');
desdir=[PathName 'power_' FileName(1:end-4)];
mkdir(desdir);
disp('Calculating power in bands...')

for j=1:length(Channels)
    Delta=zeros(length(duration),1);Alpha=Delta;Theta=Delta;
    SlowSigma=Delta;FastSigma=Delta;
    for i=1:length(duration)
        Delta(i,1)=sum(P{i,j}(find(f>=f_delta(1) & f<=f_delta(2))))*0.25;
        Alpha(i,1)=sum(P{i,j}(find(f>f_theta(1) & f<=f_theta(2))))*0.25;
        Theta(i,1)=sum(P{i,j}(find(f>f_alpha(1) & f<=f_alpha(2))))*0.25;
        SlowSigma(i,1)=sum(P{i,j}(find(f>=f_slowsigma(1) & f<=f_slowsigma(2))))*0.25;
        FastSigma(i,1)=sum(P{i,j}(find(f>=f_fastsigma(1) & f<=f_fastsigma(2))))*0.25;
    end
    output=table(Stage,duration,Delta, Alpha, Theta, SlowSigma, FastSigma);
    writetable(output,[desdir PathName(end) Channels{1,j} '.csv'],'Delimiter',',','QuoteStrings',true);
end
%     

msgbox('Done! Find the csv outputs. Checking for outliers is recommended!')
end


function [X,Channels,fs]=edf2mat_samefs(fileloc)
% [FileName,PathName] = uigetfile('/bazhlab/naji/home/EDFs_ACH_500Hz/*.edf','Select the edf data file');
% load([PathName FileName]);
fid=fopen(fileloc);% in format of [PathName FileName]
a=fread(fid,236,'*char');
ndr=fread(fid,8,'*char');
ndr=str2double(ndr'); %number of data records in sec
a=fread(fid,8,'*char');
drdur=str2double(a'); %duration of each data record in sec
ns=fread(fid,4,'*char'); ns=ns'; ns=str2double(ns);% number of signal channels
Channels=cell(1,ns);
for i=1:ns
    C=fread(fid,16,'*char');C=C';
    Channels{i}=C;
end
fread(fid,ns*80,'*char'); % channel transducer type can be extracted
fread(fid,ns*8,'*char'); %channel physical dimension can be extracted
phmn=zeros(1,ns); phmx=phmn;dmn=phmn;dmx=dmn;
for i=1:ns
    pm=fread(fid,8,'*char');pm=pm';
    phmn(i)=str2double(pm);
end                         %phys min
for i=1:ns
    pm=fread(fid,8,'*char'); pm=pm';
    phmx(i)=str2double(pm);%phys max
end
for i=1:ns
    dm=fread(fid,8,'*char');dm=dm'; 
    dmn(i)=str2double(dm);
end                         %dig min
for i=1:ns
    dx=fread(fid,8,'*char'); dx=dx';
    dmx(i)=str2double(dx);
end                         %dig max
scalefac=(phmx-phmn)./(dmx-dmn);
dc=phmx-scalefac.*dmx;

fread(fid,ns*80,'*char'); % prefilters
nr=zeros(1,ns);
for i=1:ns
    nrc=fread(fid,8,'*char'); nrc=nrc';
    nr(i)=str2double(nrc); %number of samples in each data record
end
if sum(ismember(nr,nr(1)))==length(nr)
    fs=nr(1);
else
    sprintf('Data cant be stored in a single matrix')
end
fread(fid,ns*32,'*char');
ch_fs=nr/drdur;
if mean(nr)==nr(1) && mean(ch_fs)==ch_fs(1)
X=zeros(ns,nr(1)*ndr);
% for i=1:ns
%     X{i,1}=zeros(1,nr(i)*ndr);
% end
fs=ch_fs(1);
end
disp('Reading EDF file...');
for i=1:ndr
    for j=1:ns
        s=fread(fid,nr(j),'int16').*scalefac(j)+dc(j);s=s';
        X(j,(i-1)*nr(j)+1:i*nr(j))=s;
    end

end

fclose(fid);
end

