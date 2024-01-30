function [move_trials, allspeeds, meanspeeds] = Intan_digital_movement_HW(num_trials,mmvtA,mmvtB,field_trials,runthresh,fs)
%4/19/2018 HW Add in stuff to break up fully running, fully stationary, or
%start/stop movt during trial to more specifically seperate out trials for
%2phothon

% Modified by HW 9/2017 for imec (removed built in information about
% sampling rate

%1/28/2017
%fixed issue where last speed/time point being used was taking difference
%of nonexistent zero value (use nan array instead)

% Movement decoding for the digital encoder w/ intan
% Ethan McBride - 5/5/2016
% Need to make the running/not running cutoff choice more robustly...

%field trials = on/off times per trial
% num_trials=length(trial_num);
% runthresh=0.5;
% field_trials=[on_times' off_times'];

%fs=20000; %fs=25000 imec, fs=20000 others. X20 to field trials if others
%edited to remove *20 and change fs. Don't need X20s to on/off times
%because on/off times being used are at sampling rate (not downsampled)
%%
distVect = nan(10000,3,num_trials);

for t=1:num_trials
    if rem(t,100) == 0
        fprintf('Processing trial # %d \n',t)
    end
    if field_trials(t,1)>length(mmvtA) || field_trials(t,2)>length(mmvtA)
        tempA = zeros(1,length(field_trials(t,1):field_trials(t,2)));
        tempB = zeros(1,length(field_trials(t,1):field_trials(t,2)));
        fprintf('Trial # %d has no movement data ... look into this \n',t)
    else
        tempA = mmvtA(field_trials(t,1):field_trials(t,2));
        tempB = mmvtB(field_trials(t,1):field_trials(t,2));
    end
    k=1;
    for i=2:1:length(tempA) %i= step through each sample for given on/off trial
        if tempA(i)~=tempA(i-1) %if there is a change in mmvtA
            if tempA(i)==1 %Aup
                if tempB(i)==0 %if B is down, forward motion
                    %forward
                    if k==1
                        distVect(k,1,t) = 0.25; %+ or - 0.25 because each edge is 1/4 of a blip
                    else
                        distVect(k,1,t) = distVect(k-1,1,t)+0.25;
                    end
                    distVect(k,2,t) = i/fs;
                    distVect(k,3,t) = (double(field_trials(t,1))+i)/fs;
                    k=k+1;
                elseif tempB(i)==1 %if B is up, backward motion
                    %backward
                    if k==1
                        distVect(k,1,t) = -0.25;
                    else
                        distVect(k,1,t) = distVect(k-1,1,t)-0.25;
                    end
                    distVect(k,2,t) = i/fs;
                    distVect(k,3,t) = (double(field_trials(t,1))+i)/fs;
                    k=k+1;
                end
            elseif tempA(i)==0 %Adown
                if tempB(i)==1 %if B is up, forward motion
                    %forward
                    if k==1
                        distVect(k,1,t) = 0.25;
                    else
                        distVect(k,1,t) = distVect(k-1,1,t)+0.25;
                    end
                    distVect(k,2,t) = i/fs;
                    distVect(k,3,t) = (double(field_trials(t,1))+i)/fs;
                    k=k+1;
                elseif tempB(i)==0 %if B is down, backward motion
                    %backward
                    if k==1
                        distVect(k,1,t) = -0.25;
                    else
                        distVect(k,1,t) = distVect(k-1,1,t)-0.25;
                    end
                    distVect(k,2,t) = i/fs;
                    distVect(k,3,t) = (double(field_trials(t,1))+i)/fs;
                    k=k+1;
                end
            end
        elseif tempB(i)~=tempB(i-1) %if there is a change in mmvtB
            if tempB(i)==1 %Bup
                if tempA(i)==1 %if A is up, forward motion
                    %forward
                    if k==1
                        distVect(k,1,t) = 0.25;
                    else
                        distVect(k,1,t) = distVect(k-1,1,t)+0.25;
                    end
                    distVect(k,2,t) = i/fs;
                    distVect(k,3,t) = (double(field_trials(t,1))+i)/fs;
                    k=k+1;
                elseif tempA(i)==0 %if A is down, backward motion
                    %backward
                    if k==1
                        distVect(k,1,t) = -0.25;
                    else
                        distVect(k,1,t) = distVect(k-1,1,t)-0.25;
                    end
                    distVect(k,2,t) = i/fs;
                    distVect(k,3,t) = (double(field_trials(t,1))+i)/fs;
                    k=k+1;
                end
            elseif tempB(i)==0 %Bdown
                if tempA(i)==0 %if A is down, forward motion
                    %forward
                    if k==1
                        distVect(k,1,t) = 0.25;
                    else
                        distVect(k,1,t) = distVect(k-1,1,t)+0.25;
                    end
                    distVect(k,2,t) = i/fs;
                    distVect(k,3,t) = (double(field_trials(t,1))+i)/fs;
                    k=k+1;
                elseif tempA(i)==1 %if A is up, backward motion
                    %backward
                    if k==1
                        distVect(k,1,t) = -0.25;
                    else
                        distVect(k,1,t) = distVect(k-1,1,t)-0.25;
                    end
                    distVect(k,2,t) = i/fs;
                    distVect(k,3,t) = (double(field_trials(t,1))+i)/fs;
                    k=k+1;
                end
            end
        end
    end  
end
%%
    %get true distance per up/down pulse
    distVect(:,1,:) = distVect(:,1,:)./256 .* -47.75; %if the circumference is 47.75cm %why -???
   %multiply by negative to get direction wheel moves (ie forward motion is the wheel turning backwards
   
    speed = diff(distVect(:,1,:))./diff(distVect(:,2,:)); %cm/s to get speed

    maxdistances = reshape(max(abs(distVect(:,1,:))),num_trials,1);
    maxspeeds = reshape(max(abs(speed(:,1,:))),num_trials,1);
    meanspeeds = reshape(nanmean(abs(speed(:,1,:))),num_trials,1); %abs value speed to make positive?
    
    %get time point as the time in between 2 time times (ie speed is
    %calculated from distance 2- distance 1 / (time 2- time 1), so time
    %associated with this speed is halfway in between time 1 and time 2
    speedtimepoints=reshape(distVect(2:end,3,:)-(diff(distVect(:,3,:))/2),length(speed(:,1,1)),num_trials);
    
    %allspeeds(:,:,1)= speed
    %allspeeds(:,:,2)= speed time points
    allspeeds = reshape(speed,length(speed(:,1,1)),num_trials);
    allspeeds(1:length(speed(:,1,1)),1:num_trials,2) = speedtimepoints;
    
    allspeeds(isnan(allspeeds))=0;
    meanspeeds(isnan(meanspeeds))=0;
    maxspeeds(isnan(maxspeeds))=0;
    
    move_trials = zeros(num_trials,1);
    move_trials(find(abs(meanspeeds)>=runthresh)) = 1; %run trials = any trials where run speed greater than threshold regardless of direction
    
%     figure;
%     subplot(211)
%     hist(maxdistances,20)
%     ylabel('#trials')
%     xlabel('distance traveled (cm)')
%     subplot(212)
%     hist(meanspeeds,20)
%     ylabel('#trials')
%     xlabel('mean speed (cm/sec)')