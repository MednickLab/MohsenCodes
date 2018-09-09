function out=mySOstats(edfPath,ScoresPath)

% Refer to Dang-Vu et al PNAS 2008 for more details

%                       . . b
%                     .     . 
% ----------- .<---c--><-d-->--. -------------- 
%                .    .
%                a  .

% Fixed parameters: a>80, a+b>140, 300<c<1500 ms, d<1000ms
% if you like to modify the thresholds refer to function mysodetect at the
% end of this program
% dlta: a>40, a+b>75, 300<c<1500 ms, d<1000ms

%%%%% Parameters %%%%%%%%%
q=4; % number of quartiles (quantiles)
th_min=-200; % amplitude rejection threshold for SOs
th_max=200;
epl=30; % Only change if you didn't score the sleep data in 30 sec epochs
saveout=1; % 1/0. 0 if you don't want to save output otherwise 1
%%%%%%%%%%

% load sleep scores
load(ScoresPath);
mrk=stageData.stages;
mrk(find(mrk==7))=-1;

% load edf and store in X
[X,Channels,fs]=myedf2mat_samefs(edfPath);

% find bouts and quartiles (quantiles if q<4)
t=find((mrk(2:end)-mrk(1:end-1))~=0);
smp=[0;t*epl*fs;size(mrk,1)*epl*fs]';
bnd=zeros(length(smp)-1,2);
for i=1:length(smp)-1
    bnd(i,:)=[smp(i)+1 smp(i+1)]; % beginning and end of each bout
end


o=find(mrk>0); e=find(mrk(end:-1:1)>0);

if ((length(mrk)-e(1)+1-o(1))/q)-floor((length(mrk)-e(1)+1-o(1))/q)<=0.5
    qn30=floor((length(mrk)-e(1)+1-o(1))/q);
else
    qn30=ceil((length(mrk)-e(1)+1-o(1))/q);
end

% Filtering
disp('filter_swa');
[b1,a1]=butter(4,4/(fs/2),'low');
[b2,a2]=butter(4,0.16/(fs/2),'high');
X_sw=filtfilt(b1,a1,X');
X_sw=filtfilt(b2,a2,X_sw);
X_sw=X_sw';

% SO detection

so_locs_stg2=cell(1,q);
so_locs_stg3=cell(1,q);
stg2_Qmin=zeros(1,q);
stg3_Qmin=zeros(1,q);
for ii=1:q
    disp(['Quantile #' num2str(ii)]);
    
    
    mrk_q=mrk(o(1)+(ii-1)*qn30:o(1)+ii*qn30-1,:);
    t_q=find((mrk_q(2:end)-mrk_q(1:end-1))~=0);
    smp_q=[0;t_q*epl*fs;length(mrk_q)*epl*fs]'+(o(1)+(ii-1)*qn30-1)*epl*fs;
    bnd_q=zeros(length(smp_q)-1,2);
    for i=1:length(smp_q)-1
        bnd_q(i,:)=[smp_q(i)+1 smp_q(i+1)]; % beginning and end of each bout
    end
    bnd_lbl_q=mrk_q([1;t_q+1]); % bout sleep stage
    bnd_min_q=(bnd_q(:,2)-bnd_q(:,1)+1)/(fs*60);
    
    disp('Stage 2');
    % Stage2
    i2=find(bnd_lbl_q==2);
    bnd_stg2_q=bnd_q(i2,:);
    bnd_stg2_min_q=bnd_min_q(i2);
    stg2_Qmin(1,ii)=sum(bnd_stg2_min_q);
    
    tmp2=cell(length(bnd_stg2_min_q),length(Channels));
    
    if ~isempty(bnd_stg2_min_q)
        
        for i=1:length(bnd_stg2_min_q)
            for j=1:length(Channels)
                [sw,~]=mysodetect(X_sw(j,bnd_stg2_q(i,1):bnd_stg2_q(i,2)),fs);
                sw=sw+bnd_stg2_q(i,1)*ones(size(sw))-1;
                mns=zeros(1,size(sw,1)); mnst=mns; mxs=mns;
                for s=1:size(sw,1)
                    [mm,kk]=min(X_sw(j,sw(s,1):sw(s,2)));
                    mns(s)=mm;
                    mnst(s)=sw(s,1)+kk-1;
                    mxs(s)=max(X_sw(j,sw(s,1)-fs:sw(s,2)));
                end
                
                mnst(find(mns<th_min | mxs>th_max))=[];
                tmp2{i,j}=mnst';
                %
            end
        end
    end
    so_locs_stg2{1,ii}=tmp2;
    
    
    disp('SWS');
    i3=find(bnd_lbl_q==3);
    bnd_stg3_q=bnd_q(i3,:);
    bnd_stg3_min_q=bnd_min_q(i3);
    stg3_Qmin(1,ii)=sum(bnd_stg3_min_q);
    
    tmp3=cell(length(bnd_stg3_min_q),length(Channels));
    if ~isempty(bnd_stg3_min_q)
        
        for i=1:length(bnd_stg3_min_q)
            for j=1:length(Channels)
                [sw,~]=mysodetect(X_sw(j,bnd_stg3_q(i,1):bnd_stg3_q(i,2)),fs);

                sw=sw+bnd_stg3_q(i,1)*ones(size(sw))-1;
                mns=zeros(1,size(sw,1)); mnst=mns; mxs=mns;
                for s=1:size(sw,1)
                    [mm,kk]=min(X_sw(j,sw(s,1):sw(s,2)));
                    mns(s)=mm;
                    mnst(s)=sw(s,1)+kk-1;
                    mxs(s)=max(X_sw(j,sw(s,1)-fs:sw(s,2)));
                end
                
                mnst(find(mns<th_min | mxs>th_max))=[];
                tmp3{i,j}=mnst';
            end
        end
    end
    so_locs_stg3{1,ii}=tmp3;
end

locs_stg2=cell(1,length(Channels));
locs_stg3=cell(1,length(Channels));
for j=1:length(Channels)
    locs_stg2{1,j}=[];
    for ii=1:q
        locs_stg2{1,j}=[locs_stg2{1,j};cell2mat(so_locs_stg2{1,ii}(:,j))];
        locs_stg3{1,j}=[locs_stg3{1,j};cell2mat(so_locs_stg3{1,ii}(:,j))];
    end
end
Q_soN_stg2=zeros(q,length(Channels));
Q_soN_stg3=zeros(q,length(Channels));
Q_soDns_stg2=zeros(q,length(Channels));
Q_soDns_stg3=zeros(q,length(Channels));

for ii=1:q
    if ~isempty(so_locs_stg2{1,ii})
        for j=1:length(Channels)
            Q_soN_stg2(ii,j)=length(cell2mat(so_locs_stg2{1,ii}(:,j)));
            Q_soN_stg3(ii,j)=length(cell2mat(so_locs_stg3{1,ii}(:,j)));
            Q_soDns_stg2(ii,j)=length(cell2mat(so_locs_stg2{1,ii}(:,j)))/stg2_Qmin(1,ii);
            Q_soDns_stg3(ii,j)=length(cell2mat(so_locs_stg3{1,ii}(:,j)))/stg3_Qmin(1,ii);
        end
    end
end
soN_stg2=sum(Q_soN_stg2);
soN_stg3=sum(Q_soN_stg3);
soDns_stg2=soN_stg2/sum(stg2_Qmin);
soDns_stg3=soN_stg3/sum(stg3_Qmin);

out.soN_stg2=soN_stg2;
out.soN_stg3=soN_stg3;
out.soDns_stg2=soDns_stg2;
out.soDns_stg3=soDns_stg3;
out.Q_soN_stg2=Q_soN_stg2;
out.Q_soN_stg3=Q_soN_stg3;
out.Q_soDns_stg2=Q_soDns_stg2;
out.Q_soDns_stg3=Q_soDns_stg3;
out.soTime_stg2=locs_stg2;
out.soTime_stg3=locs_stg3;
out.fs=fs;
out.channels=Channels;

if saveout==1
    sl=find(edfPath=='/' | edfPath=='\');
    mkdir(edfPath(1:sl(end)),'SO_stats');
    
    save([edfPath(1:sl(end)),'SO_stats' edfPath(sl(end)) edfPath(sl(end)+1:end-4) '_SO'],'out');
e
    soN_stg2=soN_stg2';
    soN_stg3=soN_stg3';
    soDns_stg2=soDns_stg2';
    soDns_stg3=soDns_stg3';
    Q_soN_stg2=Q_soN_stg2';
    Q_soN_stg3=Q_soN_stg3';
    Q_soDns_stg2=Q_soDns_stg2';
    Q_soDns_stg3=Q_soDns_stg3';
    T = table(soN_stg2,soN_stg3,soDns_stg2,soDns_stg3,Q_soN_stg2,Q_soN_stg3,Q_soDns_stg2,Q_soDns_stg3,...
    'RowNames',Channels');
    writetable(T,[edfPath(1:sl(end)),'SO_stats' edfPath(sl(end)) edfPath(sl(end)+1:end-4) '_SO.csv'],'WriteRowNames',true,'Delimiter',',','QuoteStrings',true);

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Functions

function [can_ind_sw,can_ind_dlta]=mysodetect(x,fs)
can_ind_sw=[];can_ind_dlta=[];
zc=[];
for i=fs+1:length(x)-fs
    if  x(i)*x(i+1)<0
        zc=[zc i+1];
    end
end
if length(zc)>1
if x(zc(1))<x(zc(1)+1)
   zc(1)=[]; % first zc should be pos to neg
end
noc=floor(length(zc)/3); %number of candidate segments
zc=zc(1:3*noc);
c1=zc(1:2:end)'; c2=zc(2:2:end)'; c3=zc(3:2:end)';
l=min([length(c1) length(c2) length(c3)]);
can_ind=[c1(1:l) c2(1:l) c3(1:l)];
cd=[can_ind(:,2)-can_ind(:,1) can_ind(:,3)-can_ind(:,2)]./fs.*1000; % sec to ms. Don't change this 1000 ! 

can_ind=can_ind(cd(:,1)>300 & cd(:,1)<1500 & cd(:,2)<1000,:);
if ~isempty(can_ind)
    [rr,~]=size(can_ind);
    ab=zeros(rr,2);
for i=1:rr
    ab(i,1)=abs(min(x(can_ind(i,1):can_ind(i,2))));
    ab(i,2)=abs(max(x(can_ind(i,2):can_ind(i,3))));
end
end
if ~isempty(can_ind)
can_ind_sw=can_ind(ab(:,1)>80 & (ab(:,1)+ab(:,2))>140,:);

can_ind_dlta=can_ind(ab(:,1)>40 & (ab(:,1)+ab(:,2))>75,:);
end
end

function [X,Channels,fs]=myedf2mat_samefs(fileloc)
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


    if i==floor(ndr/5)
        disp('Progress=20%');
    elseif i==floor(ndr*0.4)
      disp('Progress=40%');
    elseif i==floor(ndr*0.6)
      disp('Progress=60%');
    elseif i==floor(0.8*ndr)
      disp('Progress=80%');
    end
end

fclose(fid);
