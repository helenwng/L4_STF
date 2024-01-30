function [RF_center,dFoF_perstim,std_perstim,dFoF_traces,otplanes]=sbx_coarseRF(fn,fr)
%fast (~5-10 minutes) analysis to determine 2 photon RF
addpath('G:\Helen\fastRF_2PTs');
%to add:
%-get dF/F traces to sanity check that you are calculating the mean
%correctly!
%automatically plot for n planes imaged
%for 2PT S
%if nargin<1
    fn='ES23_003_001';
    fr=15.49;
%elseif nargin<2
%    fr=15.49;
%    
%end
%% load analyzer info
disp('load analyzer file: ');

[filename, pathname]=uigetfile('*.analyzer','load analyzer file: ');
[trial_nums stim_time] = looper([pathname filename]);
%%
stiminfo=readtable([fn '.txt']); % 1time, 2keyboard, 3rotB, 4rotA, 5sample, 6stim, 7PD
stiminfo=table2array(stiminfo);
figure;
for i=2:6;
    subplot(5,1,i-1); hold on;
    plot(stiminfo(:,i));
    
end
movtA=round(stiminfo(:,3)./max(stiminfo(:,3))); %normalize to binary 0,1
movtB=round(stiminfo(:,2)./max(stiminfo(:,2)));

stiminfo=stiminfo(:,5); 

%% get pulses from stiminfo
threshold=3;
pulses=find(stiminfo>threshold);  %gives you the times of where pulses start-between-end
end_edges=find(diff(pulses)>1)'; %finds indices of end edges of pulses (ie the gaps betwen pulses)
start_edges=[1 end_edges+1]; %adjacent index to diff>1 is the start of the next edge plus first start pulse
end_edges=[end_edges numel(pulses)]; %add in last end pulse that isn't included w/ difference
edges= sort([start_edges end_edges]); %first edge last edge start + end edges indices
on_off=pulses(edges); %times on/off
%%movt
alltrials=size(on_off,1)/2;
on_times=on_off(1:2:end)+1;
off_times=on_off(2:2:end)+1;
fs=1000;
[move_trials allspeeds meanspeeds]=Intan_digital_movement_HW(alltrials,movtA, movtB, [on_times off_times], 0.5,fs); %size(trial_pos,1)

%% load scanbox data
global info;
load(fn);
if info.volscan~=0
    if exist( [fn '_ot_' num2str(0,'%03.f') '.sbx'])==0
        sbxsplit(fn); %split if volumetric
%         max_idx=info.max_idx;
    end
    nslices = info.otparam(3);
    fr=fr/nslices;
else
    nslices=1;
end
nslices
%% divides up the total frames in session to which plane they match up to
frames=0:info.frame(end);
frames=frames';
frames_pp=mod(frames,nslices);

plane=nan(ceil(size(frames,1)/nslices),nslices);
for i=1:nslices
   plane(1:numel(find(frames_pp==i-1)),i)=frames(frames_pp==i-1);
end

if info.frame(1)==0 %delete first frame if at zero (not a real trial start)
    info.frame(1)=[];
    info.line(1)=[];
    info.event_id(1)=[]; 
end
on=info.frame(1:2:end);  off=info.frame(2:2:end);
disp('frame on/off/diffs'); [on off off-on];
disp(['total trials =' int2str(numel(on))]);

if min(off-on)<80
    disp('SOMETHING STILL WRONG, CHECK FRAMES');
    return;
end
%% calculate dF/F per trial for each plane and which frame per plane matches up to trial info (otplane)
splitFOV = 1;
z = squeeze(sbxread(fn,1,1));
FOVbounds = {[1:size(z,1)], [1:size(z,2)]};
% FOVbounds = {[1:size(z,1)/2], [1:size(z,2)/2]; ...
%     [1:size(z,1)/2],[(size(z,2)/2 +1) :size(z,2)]; ...
%     [(size(z,1)/2 +1) :size(z,1)],[1:size(z,2)/2]; ...
%     [(size(z,1)/2 +1) :size(z,1)],[(size(z,2)/2 +1) :size(z,2)]};

k=1; done = 0;
N=numel(on); otplanes={};

tic;
dFoF_traces=nan(ceil((off(1)-on(1))/nslices),N,nslices, size(FOVbounds,1)); %row= frames, %col= trial, %z=slice, %d4 = FOV bounds
while(~done && k<=N) %loop through each trial
    disp(['working on trial ' int2str(k)]);
    try
    for i=1:nslices
        
            if k==1
               otplanes{i}=[]; 
            end
            in=find(plane(:,i)==min(plane((plane(:,i)-on(k))>=0,i)));
            in2=find(plane(:,i)==max(plane((off(k)-plane(:,i))>=0,i)));
            otplanes{i}=vertcat(otplanes{i},[in;in2]); 
            
            if nslices>1
            fname=[fn '_ot_' num2str(i-1,'%03.f')];
            else
           fname = fn;
            end
            %z = sbxread(fname,1,1);
            q = squeeze(sbxread(fname,in-1,in2-in-1)); %-1 bc index starts at 0 
            q = uint16(smoothdata(q,3,'gaussian',4)); 

            for m = 1:size(FOVbounds,1)
       
            qds = squeeze(mean(mean(q(FOVbounds{m,1},FOVbounds{m,2},:),1),2)); %average across x and y 
            %q = squeeze(mean(mean(q(),1),2)); %average across x and y 

            dFoF_traces(1:size(qds,1),k,i,m)=(qds-mean(qds(1:round(fr*stim_time(1)))))./mean(qds(1:round(fr*stim_time(1))));
            
            meandFoF(k,i,m)=(mean(qds(round(fr*stim_time(1))+1:round(fr*(stim_time(1)+stim_time(2)))))...
                -mean(qds(1:round(fr*stim_time(1)))))/mean(qds(round(1:fr*stim_time(1))));
            end
             
            %PIXELWISE.
       %TO DO: downsample image to increase SNR of dFoF calculation?
       
           %baseline per-pixel -> too noisy
%             dFoF_FOV(k,i)=(mean(q(:,:,round(fr*stim_time(1))+1:round(fr*(stim_time(1)+stim_time(2)))),3)...
%                 -mean(q(:,:,1:round(fr*stim_time(1))),3))./mean(q(:,:,round(1:fr*stim_time(1))),3);
    
               %%base line  = avg F for entire FOV?
%              dFoF_FOV{k,i}=(mean(q(:,:,round(fr*stim_time(1))+1:round(fr*(stim_time(1)+stim_time(2)))),3)...
%                  -mean(mean(mean(q(:,:,1:round(fr*stim_time(1))),3))))./mean(mean(mean(q(:,:,1:round(fr*stim_time(1))),3)));
               dFoF_FOV{k,i}=(mean(q(:,:,round(fr*stim_time(1))+1:round(fr*(stim_time(1)+stim_time(2)))),3)...
                 -mean(mean(mean(q,3))))./mean(mean(mean(q,3)));
             
                %no baseline normalization
%               dFoF_FOV{k,i}=(mean(q(:,:,round(fr*stim_time(1))+1:round(fr*(stim_time(1)+stim_time(2)))),3)...
%                  -mean(mean(mean(q(:,:,1:round(fr*stim_time(1))),3))));
            %just average F
%               dFoF_FOV{k,i}=(mean(q(:,:,round(fr*stim_time(1))+1:round(fr*(stim_time(1)+stim_time(2)))),3));

    end
    catch
        disp(['something happened at k=' int2str(k) ' and i=' int2str(i)]); done=done+1;
    end
    k=k+1;
    toc;
end


%%
trial_nums=[trial_nums, move_trials];
%%
output=coarseRetinopathy_HW_looper;

%%
%% TO DO:
%divide FOV into grids to check that RF center does not vary across FOVs
cplot={'k','r','b','c'}
apos=unique(trial_nums);
apos(apos==500) = [];
clear apos_in
clear apos_in_r
figure(2); clf;
figure(1); clf;
for m =1:size(FOVbounds,1)
for i=1:nslices
   info.otplanes(:,i)=otplanes{i}-1; %save plane number as -1 bc sbx starts at 0 index!!
   for a=1:numel(apos)
        apos_in(:,a)= ismember(trial_nums,[apos(a),0],'rows');
        dFoF_perstim(a,:)=mean(meandFoF(apos_in(:,a),i,m));
        std_perstim(a,:)=std(meandFoF(apos_in(:,a),i,m))...
            /sqrt(length(meandFoF(apos_in(:,a),i,m)));
        
        apos_in_r(:,a)= ismember(trial_nums,[apos(a),1],'rows');
        dFoF_perstim_r(a,:)=mean(meandFoF(apos_in_r(:,a),i,m));
        std_perstim_r(a,:)=std(meandFoF(apos_in_r(:,a),i,m))...
            /sqrt(length(meandFoF(apos_in_r(:,a),i,m)));
        ntrials(a,1)=numel(find(apos_in(:,a)==1));
       ntrials(a,2)=numel(find(apos_in_r(:,a)==1));
        
        figure(2*m); subplot(2,16,a); hold on;
        slicecol = { [0.25 0.25 0.25], [1 0.5 0.5]};
        
        plot(dFoF_traces(:,apos_in(:,a),i,m),'Color',slicecol{i},'LineWidth',0.1);
   
%         plot(dFoF_traces(:,apos_in(:,a),3),'m','LineWidth',0.1);
        plot(mean(mean(dFoF_traces(:,apos_in(:,a),:,m),3),2),cplot{m},'LineWidth',2);
        xlim([0 fr*6]);  
        ylim([-0.1 0.5]);
        ycords=ylim;
        area(fr*stim_time(1):fr*(stim_time(1)+stim_time(2)), ycords(2)*ones(1, numel(fr*stim_time(1):fr*(stim_time(1)+stim_time(2)))),...
            ycords(1),'FaceColor','k','FaceAlpha',0.1,'LineStyle','none')
      
         
        figure(2*m); subplot(2,16,a+16); hold on;
        
            plot(dFoF_traces(:,apos_in_r(:,a),i,m),'Color',slicecol{i},'LineWidth',0.1);
       
        plot(mean(mean(dFoF_traces(:,apos_in_r(:,a),:,m),3),2),cplot{m},'LineWidth',2);
        xlim([0 fr*6]);
        ylim([-0.1 0.50]);
        ycords=ylim;
        area(fr*stim_time(1):fr*(stim_time(1)+stim_time(2)), ycords(2)*ones(1, numel(fr*stim_time(1):fr*(stim_time(1)+stim_time(2)))),...
            ycords(1),'FaceColor','k','FaceAlpha',0.1,'LineStyle','none')
      
   end
   
   %use either running or stat trials to find preferred horz/vert bars
   [~,horz]=max(dFoF_perstim(1:5));
   [~,vert]=max(dFoF_perstim(6:end)); vert=vert+5;
    x_pos=mean([output(vert,1),output(vert,1)]);
    y_pos=mean([output(horz,2),output(horz,2)]);
    
   [~,horzr]=max(dFoF_perstim_r(1:5));
   [~,vertr]=max(dFoF_perstim_r(6:end)); vertr=vertr+5;
    x_posr=mean([output(vertr,1),output(vertr,1)]);
    y_posr=mean([output(horzr,2),output(horzr,2)]);
    
    RF_center(i,:,1)=[horz, vert, x_pos, y_pos];
    RF_center(i,:,2)=[horzr, vertr, x_posr, y_posr];
    
    %plotting stuff
    figure(2*m-1); hold on;
    subplot(3,nslices,i); hold on;
    errorbar(dFoF_perstim,std_perstim,'ko');
    errorbar(dFoF_perstim_r,std_perstim_r,'ro');
    title('dF/F per bars');
    if i==nslices
    legend('stationary','running','Location','SouthEast');
    end
    
    subplot(3,nslices,i+nslices); 
    fill([0 0 1280 1280], ... 
        [output(horz,2)-output(horz,4)/2 output(horz,2)+output(horz,4)/2 ...
        output(horz,2)+output(horz,4)/2 output(horz,2)-output(horz,4)/2] ,'b'...
        ,'FaceAlpha',0.4);
    
    hold on;
    
    %stat
    fill([output(vert,1)-output(vert,3)/2 output(vert,1)+output(vert,3)/2 ...
        output(vert,1)+output(vert,3)/2 output(vert,1)-output(vert,3)/2], ...
        [0 0 720 720],'r','FaceAlpha',0.4);
    plot(x_pos,y_pos,'kx','MarkerSize',20);
    text(x_pos+20,y_pos-15,['(' int2str(x_pos) ',' int2str(y_pos) ')']);
    xlim([0 1280]);
    ylim([0 720]);
    title(['RF for plane stationary=' int2str(i)]);
    set(gca, 'Ydir', 'reverse')
    
    %running
    subplot(3,nslices,i+2*nslices); 
    fill([0 0 1280 1280], ... 
        [output(horzr,2)-output(horzr,4)/2 output(horzr,2)+output(horzr,4)/2 ...
        output(horzr,2)+output(horzr,4)/2 output(horzr,2)-output(horzr,4)/2] ,'y'...
        ,'FaceAlpha',0.4);
    
    hold on;
    
    fill([output(vertr,1)-output(vertr,3)/2 output(vertr,1)+output(vertr,3)/2 ...
        output(vertr,1)+output(vertr,3)/2 output(vertr,1)-output(vertr,3)/2], ...
        [0 0 720 720],'g','FaceAlpha',0.4);
    plot(x_posr,y_posr,'kx','MarkerSize',20);
    text(x_posr+20,y_posr-15,['(' int2str(x_posr) ',' int2str(y_posr) ')']);
    xlim([0 1280]);
    ylim([0 720]);
    title(['RF for plane running =' int2str(i)]);
    set(gca, 'Ydir', 'reverse')

end
 
save([fn '_info.mat'],'info');
end
%%
%%%
%% compute map for horizontal + vertical preferences.
%color cold a vector -> preferred position + magnitude?
%average each pixel by weight of each position?
%calculate dF/F per position then pseudo color a color = position, and
%overlay?
tic
for a=1:numel(apos)
    apos_in(:,a)= ismember(trial_nums(:,1),[apos(a)],'rows');

   % stimavg(:,:,a) = mean(cell2mat(reshape(dFoF_FOV(trial_nums(:,1) ==apos(a)),1,1,[])),3);
   stimavg(:,:,a) = imgaussfilt(mean(cell2mat(reshape(dFoF_FOV(apos_in(:,a),1),1,1,[])),3),1);
end
stimavg= stimavg./max(max(max(stimavg)));

toc

%%
figure; hold on;

y_ps= jet(6);
for a=1:6
    ax(a) = subplot(1,6,a); hold on;
% if a==1
%         ax(1) = subplot(1,1,1); hold on;
%     else
%         ax(a) =axes('Position', ax(1).Position);
% end
%imagesc(ones(size(stimavg(:,:,1))),'AlphaData',stimavg(:,:,a));
imagesc(stimavg(:,:,a));
caxis([0 0.5]);
c = gray;
c = c.*y_ps(a,:);
colormap(ax(a),c);
set(ax(a),'ydir','reverse');
end
%%
%% rainbow overlap w/ alpha value for each color == intensity
figure; hold on;

y_ps= jet(6);
for a=1:6
    %ax(a) = subplot(1,6,a); hold on;
if a==1 
        ax(1) = subplot(1,1,1); hold on;
    else
        ax(a) =axes('Position', ax(1).Position);
end
imagesc(ones(size(stimavg(:,:,1))),'AlphaData',stimavg(:,:,a));
%imagesc(stimavg(:,:,a));
if a~=1 
set(ax(a),'color','none','visible','off');
end
set(ax(a),'ydir','reverse');
caxis([0 0.5]);
c = gray;
c = c.*y_ps(a,:);
colormap(ax(a),c);
end
linkaxes([ax(1), ax(2), ax(3), ax(4), ax(5), ax(6)]);
%%
figure; hold on;
y_ps= jet(6);
for a=1:6

        ax(a) = subplot(1,6,a); hold on;

%imagesc(ones(size(stimavg(:,:,1))),'AlphaData',stimavg(:,:,a));
imagesc(stimavg(:,:,a));
caxis([0 1]);
c = gray;
set(ax(a),'ydir','reverse');
colormap(ax(a),c);
end
%%
horz = sort(stimavg(:,:,7:16),3);
[mag,val]=max(stimavg(:,:,7:16),[],3); %preferred stim
%need a "baseline" for bottom of tuning curve
bl =mean(horz(:,:,1:4),3); %min(stimavg(:,:,7:16),[],3);%"baseline"
resp_mag = (mag-bl)./bl;
figure; 
ax1= subplot(1,2,1);
imagesc(ax1,mean(uint16(smoothdata(squeeze(sbxread(fname,0,100)),3,'gaussian',1)),3));


caxis([0 1000000]);%-1 bc index starts at 0 
ax2 =axes('Position', ax1.Position);

%convolve image to smooth it out?
%K = (1/36)*ones(6);
%imagesc(ax2,val,'AlphaData',conv2(resp_mag./prctile(resp_mag(:),90),K,'same'));

%imagesc(ax2,val,'AlphaData',mag./max(mag(:)));clear all;

imagesc(ax2,val,'AlphaData',mag*2);
%link axes 
linkaxes([ax1,ax2]) 
%%Hide the top axes 
ax2.Visible = 'off'; 
ax2.XTick = []; 
ax2.YTick = []; 
%add differenct colormap to different data if you wish 
colormap(ax1,'gray') 
colormap(ax2,'jet') 
colorbar(ax2,'Position', [0.48, 0.11, 0.01, 0.82]);
title(ax1,'horizontal RF (x)');


[mag,val]=max(stimavg(:,:,1:6),[],3);

ax3=subplot(1,2,2);
imagesc(ax3,mean(uint16(smoothdata(squeeze(sbxread(fname,0,100)),3,'gaussian',1)),3));
colormap('gray');
caxis([0 1000000]);%-1 bc index starts at 0 
ax4 =axes('Position', ax3.Position);
imagesc(ax4,val,'AlphaData',mag);
colormap('jet')
%link axes 
linkaxes([ax3,ax4]) 
%%Hide the top axes 
ax4.Visible = 'off'; 
ax4.XTick = []; 
ax4.YTick = []; 
%add differenct colormap to different data if you wish 
colormap(ax3,'gray') 
colormap(ax4,'jet') 
colorbar(ax4,'Position', [0.92, 0.11, 0.01, 0.82],'XTick',[1:1:6]);

linkaxes([ax1,ax2,ax3,ax4]);
title(ax3,'vert RF (y)');

set(gcf,'Position',[100 400 1500 500]);

%%

