function [trial_num stim_time] = looper(fa,stf_stim)
%Looper file

if nargin<2
    stf_stim=1;
end
load (fa, '-mat');

bflag = 0; %are there blank trials?
if strcmp(Analyzer.loops.conds{length(Analyzer.loops.conds)}.symbol{1},'blank')
    bflag = 1;
    nb=1;
    lb=0;
end

if any(strcmp(Analyzer.loops.conds{length(Analyzer.loops.conds)}.symbol,'light_bit'))
    lb = 1;
    nb=2;
end

cond=length(Analyzer.loops.conds); %total number of combinations of conditions (last cond is blank condition)
if bflag ==1 %if so, here's how you deal with it
    for i=1:length(Analyzer.loops.conds)-nb %# of blank trial types (ie lb on/off)
        if i==1
            if lb==1 %if there is a lightbit
                trial_num=500*ones(length(Analyzer.loops.conds{i}.symbol),...
                    (length(Analyzer.loops.conds)-2)*length(Analyzer.loops.conds{i}.repeats)...
                    +(length(Analyzer.loops.conds{cond-1}.repeats))+(length(Analyzer.loops.conds{cond}.repeats)))'; 
            else
                trial_num=500*ones(length(Analyzer.loops.conds{i}.symbol),(length(Analyzer.loops.conds)-1)*length(Analyzer.loops.conds{i}.repeats)+(length(Analyzer.loops.conds{cond}.repeats)))'; %9x20=180 even tho 175 total, bc 20 repeats for stim + 15 blanks
            end
         end 
          for k=1:length(Analyzer.loops.conds{i}.symbol)
              trial_vals(:,k)= Analyzer.loops.conds{i}.val{k};
          end
          for j=1:length(Analyzer.loops.conds{i}.repeats)
              aux_trial=Analyzer.loops.conds{i}.repeats{j}.trialno;
              trial_num(aux_trial,:)=trial_vals;
          end 
    end
    if lb==1 %go back and add in 
    for i=cond-1:cond
        if isempty(Analyzer.loops.conds{i}.val{end})
        bb=0;
        else
        bb=Analyzer.loops.conds{i}.val{end};
        end
        for j=1:length(Analyzer.loops.conds{i}.repeats)
          aux_trial=Analyzer.loops.conds{i}.repeats{j}.trialno;
          trial_num(aux_trial,end)=bb;
        end     
    end       
    end
else %if there's no blank trials, use all conds:
for i=1:length(Analyzer.loops.conds)
   if i==1;trial_num=ones(length(Analyzer.loops.conds{i}.symbol),length(Analyzer.loops.conds)*length(Analyzer.loops.conds{i}.repeats))';end ;
      for k=1:length(Analyzer.loops.conds{i}.symbol)
          trial_vals(:,k)= Analyzer.loops.conds{i}.val{k};
      end
      for j=1:length(Analyzer.loops.conds{i}.repeats)
          aux_trial=Analyzer.loops.conds{i}.repeats{j}.trialno;
          trial_num(aux_trial,:)=trial_vals;
      end 
end
end

predelay = cell2mat(Analyzer.P.param{1}(3));
postdelay = cell2mat(Analyzer.P.param{2}(3));
stim_time = cell2mat(Analyzer.P.param{3}(3));
stim_time = [predelay stim_time postdelay];


%if analyzer file is from a sf, tf, ori presenting stimulus set, find the
%blanks and set to 500 bc 1 is too similar to a real value.
%[r,c]=find(trial_num==1); %this takes all blank trials (if there were any) and sets it to 256

% if stf_stim==1
%     trial_num(r,:)=500; % sort either symbol by 256 for blank trials 
% end
end
