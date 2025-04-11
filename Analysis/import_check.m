% import check - this version ADVICE CONTEXT 

% cd and import
cd('/Users/164376/Documents/GitHub/js_ArtificialDotsAdvice/Analysis/');
cd('Raw_data')

%% pseudo code.
% load csv as .mat
% delete unused columns. 


mytab=readtable("experiment-data (52).csv");
%% START DATA WRANGLING: 

% delete the bloat columns, single entries we dont need. 
del_columns={'success',...
    'timeout',...
    'failed_images',...
    'failed_video',...
    'failed_audio',...
    'internal_node_id'};
for icol= 1:length(del_columns);

try mytab.([del_columns{icol}])=[];
catch; disp(['no "' del_columns{icol} '" column']);
end
end
%% convert to binary most true-false columns: 
forcols ={'correct', 'correct_1', 'adviceCorrect', 'correct_2'};
for icol=1:length(forcols);

    % extract the data.
    tmp = nan(height(mytab),1);
    for itrial = 1:height(mytab);
        tmp= mytab.([forcols{icol}])(itrial);
         if contains(tmp, 'true');
            tmpdata(itrial,1)= 1;
         elseif contains(tmp, 'false');
            tmpdata(itrial,1)= 0;
         else
             tmpdata(itrial,1)= nan;
         end
    end
    % rename and save after relevant data.:
    mytab.([forcols{icol} '_binary'])= tmpdata;
    mytab= movevars(mytab, [forcols{icol} '_binary'], 'after', forcols{icol});

end



%% there are various instructions and practice screens/trials recorded in the output. 
% the experiment proper begins at the first 'stimulus==block message'

%first sanity check, staircase performance during practice

exp_startrow = find(contains(mytab.stimulus,'expInstruc'), 1,'last');
% only one trial type in practice (confidence_first)
conf1_rows= find(contains(mytab.task(1:exp_startrow), 'confidence_first'));
practab = mytab(conf1_rows,:);


%% now plot for sanity check:
clf;
ntrials = 1:height(practab);
subplot(221);
plot(ntrials, practab.threshold);
ylabel('dotdifference per trial');
shg

% rolling average
estAv = nan(length(ntrials),1);
for itrial = ntrials
estAv(itrial,1)= sum(practab.correct_binary(1:itrial))/ itrial;
end
%%
subplot(222);
plot(ntrials, estAv, 'r-o')
ylabel('rolling accuracy');
shg
%%
%% now for the experiment.
consolidated_rows = find(contains(mytab.stimulus,'consolidated-trial-data'));
ctable = mytab(consolidated_rows,:);
%%
% need to do some reordering. 
% delete the bloat columns, single entries we dont need. 
del_columns={'trial_type',...
    'trial_index',...
    'participant_age',...
    'participant_gender',...
    'rt',...
    'stimulus',...
    'response',...
    'task',...
    'pixels_per_unit',...
    'itrial','iblock', 'leftDots', 'rightDots','difference',...
    'selected_side','confidence_value', 'left_confidence', 'right_confidence',...
    'correct', 'threshold', 'block_num', 'advice_direction', 'top_prompt','bottom_prompt',...
    'top_response','bottom_response'};


for icol= 1:length(del_columns);

try ctable.([del_columns{icol}])=[];
catch; disp(['no "' del_columns{icol} '" column']);
end
end
%tidy the order. 
ctable= movevars(ctable,'block','After',1);
ctable= movevars(ctable,'trial','After',2);

%% seems to be an indexing bug in the consolidated data rows. adviceConds 
% are 1,2,3 (choice, forced, no advice). Some 2 and 3 are displaying the
% choice selected from the previous '1'. 
% change adviceChoice to nan for adviceConds 2,3.

for itrial= 1:height(ctable)
if ctable.adviceCond(itrial)>1
    ctable.adviceChoice{itrial}='NaN'; % text to make string search easier
elseif ctable.adviceCond(itrial)==1
    % do nothing for now.

end
end

%% include a column for whether advice was congruent with decision 1.
adviceCongr=[];
for itrial= 1:height(ctable);
dir1 = ctable.selectedSide_first(itrial);

if contains(ctable.adviceDir(itrial),dir1);
    adviceCongr(itrial,1)=1;
else
    adviceCongr(itrial,1)=0;
end
end
%add to table and position appropr.
ctable.adviceCongruent = adviceCongr;
ctable= movevars(ctable, 'adviceCongruent', 'After', 'adviceDir');

%%  include integer (1,2) for left and right to make data processing easier.
tmp=nan;
lr= find(contains(ctable.selectedSide_first,'left'));
rr= find(contains(ctable.selectedSide_first,'right'));
tmp(lr,1) = 1;
tmp(rr,1) = 2;
ctable.selectedSide_int_first = tmp;
ctable= movevars(ctable, 'selectedSide_int_first', 'After', 'selectedSide_first');

% same for second:
tmp=nan;
lr= find(contains(ctable.selectedSide_second,'left'));
rr= find(contains(ctable.selectedSide_second,'right'));
tmp(lr,1) = 1;
tmp(rr,1) = 2;
ctable.selectedSide_int_second = tmp;
ctable= movevars(ctable, 'selectedSide_int_second', 'After', 'selectedSide_second');
%% now we can confirm basic stats:
% staircase first (total exp). convert correct_1 to binary:
% convert to binary:
[tmpdata, tmpdata2]= deal(nan(height(ctable),1));

for itrial = 1:height(ctable);

    tmp = ctable.correct_1(itrial);
    if contains(tmp, 'true');
        tmpdata(itrial,1)= 1;
    else
        tmpdata(itrial,1)= 0;
    end
%
    tmp = ctable.correct_2(itrial);
    if contains(tmp, 'true');
        tmpdata2(itrial,1)= 1;
    else
        tmpdata2(itrial,1)= 0;
    end
end
ctable.correct_1binary = tmpdata;
ctable.correct_2binary = tmpdata2;
ctable= movevars(ctable,'correct_1binary','After','correct_1');
ctable= movevars(ctable,'correct_2binary','After','correct_2');
%% plot total exp accuracy and thresholds.

clf;
ntrials = 1:height(ctable);
subplot(221);
plot(ntrials, ctable.dotsDiff);
ylabel('dotdifference per trial');
shg

% rolling average
[estAv, estAv2] = deal(nan(length(ntrials),1));

for itrial = ntrials
estAv(itrial,1)= sum(ctable.correct_1binary(1:itrial))/ itrial;
estAv2(itrial,1)= sum(ctable.correct_2binary(1:itrial))/ itrial;
end

%
subplot(222);
plot(ntrials, estAv, 'b-o'); hold on;
plot(ntrials, estAv2, 'r-o'); hold on;
ylabel('rolling accuracy');
legend('first resp', 'second');
shg

%% what is the total accuracy on c1, then c2 by option? 
accC1 = sum(ctable.correct_1binary)/ height(ctable);
accC2 = sum(ctable.correct_2binary)/ height(ctable);


choiceY_idx = intersect(find(ctable.adviceCond==1), find(contains(ctable.adviceChoice,'yes')));
choiceN_idx = intersect(find(ctable.adviceCond==1), find(contains(ctable.adviceChoice,'no')));
forced_idx = find(ctable.adviceCond==2);
noAdv_idx = find(ctable.adviceCond==3);

subplot(2,2,3);
bar(1:2, [accC1, accC2]);
ylabel('Accuract on first and second choice');
ylabel('proportion correct');
set(gca,'xticklabels', {'first', 'second'});
xlabel('response order');
ylim([.5 1]);

%acc per type
chY_acc = sum(ctable.correct_2binary(choiceY_idx))/ length(choiceY_idx);
chN_acc = sum(ctable.correct_2binary(choiceN_idx))/ length(choiceN_idx);
forced_acc= sum(ctable.correct_2binary(forced_idx))/ length(forced_idx);
noAdv_acc= sum(ctable.correct_2binary(noAdv_idx))/ length(noAdv_idx);

subplot(224);
bar(1:4, [chY_acc, chN_acc, forced_acc, noAdv_acc]);
xlabel('advice condition');
ylabel('proportion correct');
set(gca,'xticklabels', {'chose Yes', 'chose No', 'forced', 'no advice'})
shg


%% does confidence track accuracy?
respwas ={'first', 'second'};
PFX_Conf=[];
PFX_RT=[];
for iresp = 1:2
    corr_idx = find(ctable.(['correct_' num2str(iresp) 'binary'])==1);
    err_idx = find(ctable.(['correct_' num2str(iresp) 'binary'])==0);

    % collect responses.
    confData = ctable.(['confidence_' respwas{iresp}]);
    rtData= ctable.(['rt_' num2str(iresp)]);

    meanConf = [nanmean(confData(corr_idx)),...
        nanmean(confData(err_idx))];

    meanRT = [nanmean(rtData(corr_idx)), ...
        nanmean(rtData(err_idx))];

    PFX_Conf(iresp,1:2) = meanConf;
    PFX_RT(iresp,1:2) = meanRT;
end
%%plot

subplot(2,2,3);
bar(PFX_Conf); ylabel('mean confidence')
set(gca,'xticklabels', {'correct', 'error'});
legend('first conf', 'second conf', 'Location', 'northoutside')

subplot(2,2,4); 
bar(PFX_RT);ylabel('mean RTs');
% ylim([0 1]);
set(gca,'xticklabels', {'correct', 'error'});
legend('first conf', 'second conf', 'Location', 'northoutside')
shg
%% next figure. plot agreement:
% pseudo:
% Proportion selected side = advice side (?)
% Confidence change (congr vs incongr)

%prop selected
choiceTrials= find(ctable.adviceCond==1);
choices = ctable.adviceChoice(choiceTrials);
choseY = find(contains(choices, 'yes'));
choseN = find(contains(choices, 'no'));
propY= length(choseY)/length(choiceTrials);
propN= length(choseN)/length(choiceTrials);
% initial Confidence for each. 
conf1onchoice = ctable.confidence_first(choiceTrials)
conf_chY = conf1onchoice(choseY);
conf_chN = conf1onchoice(choseN);
clf;
subplot(1,3,1)
bar(1, [propY;propN], 'stacked');
legend('yes', 'no');
ylabel('proportion'); title('would you like advice?')

subplot(132); 
bar(1:2, [mean(conf_chY), mean(conf_chN)]); ylim([50 100]);
title('initial confidence'); set(gca,'xticklabels', {'yes', 'no'}); xlabel('Would you like advice?');
ylabel('initial confidence')

% now final confidence:
% initial Confidence for each. 
conf2onchoice = ctable.confidence_second(choiceTrials);
conf_chY = conf2onchoice(choseY);
conf_chN = conf2onchoice(choseN);

subplot(133); 
bar(1:2, [mean(conf_chY), mean(conf_chN)]); ylim([50 100]);
title('final confidence'); set(gca,'xticklabels', {'yes', 'no'}); xlabel('Would you like advice?');
ylabel('final confidence')
%% accuracy when advice congr/not
congrtrials = find(ctable.adviceCongruent==1);
incongrtrials = find(ctable.adviceCongruent==0);



accData = [nansum(ctable.correct_2binary(congrtrials))./length(congrtrials),...
            nansum(ctable.correct_2binary(incongrtrials))./length(incongrtrials)];
% confidence
%%

subplot(2,2,1);
bar(bothDSB); ylabel('mean DSB (complete)'); 
set(gca,'xticklabels', {'preDSB', 'postDSB'});
legend('easy', 'hard');
ylim([0 1])

subplot(2,2,4);
bar(bothDSB_digits); ylabel('mean DSB (by digits)'); 
set(gca,'xticklabels', {'preDSB', 'postDSB'});
legend('easy', 'hard');
ylim([0 1])

% subplot(234);
% bar(postDSB); ylabel('mean postDSB')
shg
% plot()

%% now to calculate the pre and post by image type.
% (vegetation level)

barData=[];
v_1 = find(ctable.vegetationLevel==1);
v_2 = find(ctable.vegetationLevel==2);
v_3 = find(ctable.vegetationLevel==3);


clf;
meanRT = [nanmean(ctable.Vis_meanRT(v_1)),nanmean(ctable.Vis_meanRT(v_2)), nanmean(ctable.Vis_meanRT(v_3))];
meanAcc = [nanmean(ctable.Vis_propCorr(v_1)),nanmean(ctable.Vis_propCorr(v_2)),nanmean(ctable.Vis_propCorr(v_3))];

preDSB = [nanmean(ctable.DSBpre_propcorrect(v_1)),nanmean(ctable.DSBpre_propcorrect(v_2)),nanmean(ctable.DSBpre_propcorrect(v_3))];
preDSB_digits = [nanmean(ctable.DSBpre_digitaccuracy(v_1)), nanmean(ctable.DSBpre_digitaccuracy(v_2)),nanmean(ctable.DSBpre_digitaccuracy(v_3))];

postDSB= [nanmean(ctable.DSBpost_propcorrect(v_1)),nanmean(ctable.DSBpost_propcorrect(v_2)), nanmean(ctable.DSBpost_propcorrect(v_3))];

postDSB_digits = [nanmean(ctable.DSBpost_digitaccuracy(v_1)),nanmean(ctable.DSBpost_digitaccuracy(v_2)),nanmean(ctable.DSBpost_digitaccuracy(v_3))];

bothDSB= [preDSB;postDSB];
bothDSB_digits= [preDSB_digits; postDSB_digits];

subplot(2,2,1);
bar(meanRT); ylabel('mean vissearch RT')
set(gca,'xticklabels', {'veg1', 'veg2', 'veg3'});

subplot(2,2,2); 
bar(meanAcc);ylabel('mean vissearch Acc');
ylim([0 1]);
set(gca,'xticklabels', {'veg1', 'veg2', 'veg3'});

subplot(2,2,3);
bar(bothDSB); ylabel('mean DSB (complete)'); 
set(gca,'xticklabels', {'preDSB', 'postDSB'});
legend('veg1', 'veg2', 'veg3');
ylim([0 1])

subplot(2,2,4);
bar(bothDSB_digits); ylabel('mean DSB (by digits)'); 
set(gca,'xticklabels', {'preDSB', 'postDSB'});
legend('veg1', 'veg2', 'veg3');
ylim([0 1])
shg