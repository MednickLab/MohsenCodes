% Developed by Mohsen Naji on Aug 4, 2017 for Sleep and Cognition lab, UCR
% mohsen_dot_bme_at_gmail_dot_com
% I assume the data is re-referenced
% Ignore the outputs for E1/E2, CHIN/EMG and ECG channels
function myeegpwr_folder_v2()
clearvars
% ctu=0;
choice = questdlg('Are you a Windows or Mac/linux user?', 'Select one','Windows', 'Mac/Linux','Mac/Linux');
if choice(1)=='M'
    sb='/';
else
    sb='\';
end

disp('1) select folder for edf files; 2) select folder for marker files');
disp('This program assumes edf files contain re-referenced data(to mastoids or earlobes)');

foldernameE = uigetdir(pwd,'edf folder');
allfilesE = dir(fullfile(foldernameE));
is=find([allfilesE.isdir]==0);
foldernameM = uigetdir(pwd,'marker folder');

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

for sj=1:length(is)
    name=allfilesE(is(sj)).name;
if name(1)~='.'
        
        load([foldernameM sb name(1:end-3)]);
        
    disp(['Marker is loaded for ' name(1:end-3)]);
    
[X,Channels,fs]=edf2mat_samefs([foldernameE sb name]);
disp(['EDF data is loaded for ' name(1:end-3)]);

wins=stageData.win;
mrk=stageData.stages;
t=find((mrk(2:end)-mrk(1:end-1))~=0);
smp=[0;t*30*fs;size(mrk,1)*30*fs]';
bnd=zeros(length(smp)-1,2);
for i=1:length(smp)-1
    bnd(i,:)=[smp(i)+1 smp(i+1)]; % beginning and end of each bout
end
Stage=mrk([1;t+1]); % bout sleep stage
duration=(bnd(:,2)-bnd(:,1)+1)/(fs*60); % bout duration in minute
P=cell(length(duration)-1,length(Channels));
spins={'\\','|','/','-'};
reverseStr = 'Power spectrum estimation ';
for i=1:length(duration)-1
    si=mod(i,4); 
    if si==0
        si=4;
    end
    msg = (spins{si});
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, 1);
    
    x=X(:,bnd(i,1):bnd(i,2));
    [pxx,f]=pwelch(detrend(x'),hanning(4*fs),2*fs,0:1/4:(fs/2),fs);
    for j=1:length(Channels)
        P{i,j}=pxx(:,j);
    end
end
Stage(end)=[];
duration(end)=[];
fprintf('\n');
desdir=[foldernameE sb 'power_' name(1:end-4)];
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
    writetable(output,[desdir sb Channels{1,j} '.csv'],'Delimiter',',','QuoteStrings',true);
    
    wake_delta=Delta(find(Stage==0)); P_wake_delta(j,1)=mean(wake_delta(find(wake_delta<(mean(wake_delta)+3*std(wake_delta)))));
    stg1_delta=Delta(find(Stage==1)); P_stg1_delta(j,1)=mean(stg1_delta(find(stg1_delta<(mean(stg1_delta)+3*std(stg1_delta)))));
    stg2_delta=Delta(find(Stage==2)); P_stg2_delta(j,1)=mean(stg2_delta(find(stg2_delta<(mean(stg2_delta)+3*std(stg2_delta)))));
    stg3_delta=Delta(find(Stage==3)); P_stg3_delta(j,1)=mean(stg3_delta(find(stg3_delta<(mean(stg3_delta)+3*std(stg3_delta)))));
    rem_delta=Delta(find(Stage==5)); P_rem_delta(j,1)=mean(rem_delta(find(rem_delta<(mean(rem_delta)+3*std(rem_delta)))));
    
    wake_theta=Theta(find(Stage==0)); P_wake_theta(j,1)=mean(wake_theta(find(wake_theta<(mean(wake_theta)+3*std(wake_theta)))));
    stg1_theta=Theta(find(Stage==1)); P_stg1_theta(j,1)=mean(stg1_theta(find(stg1_theta<(mean(stg1_theta)+3*std(stg1_theta)))));
    stg2_theta=Theta(find(Stage==2)); P_stg2_theta(j,1)=mean(stg2_theta(find(stg2_theta<(mean(stg2_theta)+3*std(stg2_theta)))));
    stg3_theta=Theta(find(Stage==3)); P_stg3_theta(j,1)=mean(stg3_theta(find(stg3_theta<(mean(stg3_theta)+3*std(stg3_theta)))));
    rem_theta=Theta(find(Stage==5)); P_rem_theta(j,1)=mean(rem_theta(find(rem_theta<(mean(rem_theta)+3*std(rem_theta)))));
    
    wake_alpha=Alpha(find(Stage==0)); P_wake_alpha(j,1)=mean(wake_alpha(find(wake_alpha<(mean(wake_alpha)+3*std(wake_alpha)))));
    stg1_alpha=Alpha(find(Stage==1)); P_stg1_alpha(j,1)=mean(stg1_alpha(find(stg1_alpha<(mean(stg1_alpha)+3*std(stg1_alpha)))));
    stg2_alpha=Alpha(find(Stage==2)); P_stg2_alpha(j,1)=mean(stg2_alpha(find(stg2_alpha<(mean(stg2_alpha)+3*std(stg2_alpha)))));
    stg3_alpha=Alpha(find(Stage==3)); P_stg3_alpha(j,1)=mean(stg3_alpha(find(stg3_alpha<(mean(stg3_alpha)+3*std(stg3_alpha)))));
    rem_alpha=Alpha(find(Stage==5)); P_rem_alpha(j,1)=mean(rem_alpha(find(rem_alpha<(mean(rem_alpha)+3*std(rem_alpha)))));
    
    wake_slowsigma=SlowSigma(find(Stage==0)); P_wake_slowsigma(j,1)=mean(wake_slowsigma(find(wake_slowsigma<(mean(wake_slowsigma)+3*std(wake_slowsigma)))));
    stg1_slowsigma=SlowSigma(find(Stage==1)); P_stg1_slowsigma(j,1)=mean(stg1_slowsigma(find(stg1_slowsigma<(mean(stg1_slowsigma)+3*std(stg1_slowsigma)))));
    stg2_slowsigma=SlowSigma(find(Stage==2)); P_stg2_slowsigma(j,1)=mean(stg2_slowsigma(find(stg2_slowsigma<(mean(stg2_slowsigma)+3*std(stg2_slowsigma)))));
    stg3_slowsigma=SlowSigma(find(Stage==3)); P_stg3_slowsigma(j,1)=mean(stg3_slowsigma(find(stg3_slowsigma<(mean(stg3_slowsigma)+3*std(stg3_slowsigma)))));
    rem_slowsigma=SlowSigma(find(Stage==5)); P_rem_slowsigma(j,1)=mean(rem_slowsigma(find(rem_slowsigma<(mean(rem_slowsigma)+3*std(rem_slowsigma)))));
    
    wake_fastsigma=FastSigma(find(Stage==0)); P_wake_fastsigma(j,1)=mean(wake_fastsigma(find(wake_fastsigma<(mean(wake_fastsigma)+3*std(wake_fastsigma)))));
    stg1_fastsigma=FastSigma(find(Stage==1)); P_stg1_fastsigma(j,1)=mean(stg1_fastsigma(find(stg1_fastsigma<(mean(stg1_fastsigma)+3*std(stg1_fastsigma)))));
    stg2_fastsigma=FastSigma(find(Stage==2)); P_stg2_fastsigma(j,1)=mean(stg2_fastsigma(find(stg2_fastsigma<(mean(stg2_fastsigma)+3*std(stg2_fastsigma)))));
    stg3_fastsigma=FastSigma(find(Stage==3)); P_stg3_fastsigma(j,1)=mean(stg3_fastsigma(find(stg3_fastsigma<(mean(stg3_fastsigma)+3*std(stg3_fastsigma)))));
    rem_fastsigma=FastSigma(find(Stage==5)); P_rem_fastsigma(j,1)=mean(rem_fastsigma(find(rem_fastsigma<(mean(rem_fastsigma)+3*std(rem_fastsigma)))));
    
end
outputA=table(Channels',P_wake_delta,P_wake_theta,P_wake_alpha,P_wake_slowsigma,P_wake_fastsigma,...
    P_stg1_delta,P_stg1_theta,P_stg1_alpha,P_stg1_slowsigma,P_stg1_fastsigma,...
    P_stg2_delta,P_stg2_theta,P_stg2_alpha,P_stg2_slowsigma,P_stg2_fastsigma,...
    P_stg3_delta,P_stg3_theta,P_stg3_alpha,P_stg3_slowsigma,P_stg3_fastsigma,...
    P_rem_delta,P_rem_theta,P_rem_alpha,P_rem_slowsigma,P_rem_fastsigma);
    writetable(outputA,[desdir sb 'averages.csv'],'Delimiter',',','QuoteStrings',true);
end
end
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

