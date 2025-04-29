% import check - this version ADVICE CONTEXT


clear all;
close all;
jobs=[];

jobs.import_andwrangle = 1; % resaves as matlab tables.

jobs.plot_PP_summary=1;
jobs.createGFX_table=1;
nCreate = 10; 

jobs.plotGFX=1;

  % cd and import
  cd('/Users/164376/Documents/GitHub/js_ArtificialDotsAdvice/Analysis/');
  cd('Raw_data')
  datadir=pwd;

pfols = dir([pwd filesep 'experiment-data*']);
disp([num2str(length(pfols)) ', datafiles detected'])
%%
if jobs.import_andwrangle==1

%%
    for ippant=1:length(pfols)
        % cd and import
        cd('/Users/164376/Documents/GitHub/js_ArtificialDotsAdvice/Analysis/');
        cd('Raw_data')

        if ippant <10
            pnum= ['0' num2str(ippant)]; % maintains order 01,02,03...,10..
        else
            pnum = [num2str(ippant)];
        end
        % pseudo code.
        % load csv as .mat
        % delete unused columns.


        % mytab=readtable("experiment-data-MD.csv");
        mytab=readtable(pfols(ippant).name);


        %%  START DATA WRANGLING:
        % 1) delete unnecessary columns
        % 2) correct for the bug (all nans for correct_2 in 1 dataset).
        % 3) convert to binary most true-false columns (for easy accuracy
        % summaries).
        % 4) convert all conf from -100 -50, 50 100 to -50 0 0 50, then zscored.
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

%% %% seems the conversion to matlab struggles when parsing the JSON strings that
% make the responses to survey data (PRS etc). The data is visible in the
% CSV, but not mytab.

% Read the CSV file using fopen/textscan for more control
fid = fopen(pfols(ippant).name);

% Read header line to get column names
headerLine = fgetl(fid);
headers = textscan(headerLine, '%s', 'Delimiter', ',');
headers = headers{1};

% Find the index of the important columns
responseColIdx = find(contains(headers, 'response'));
taskColIdx = find(contains(headers, 'task'));
attnCheckIdx = find(contains(headers, 'attention'));
% Read the rest of the file as text
data = textscan(fid, '%s', 'Delimiter', '\n');
data = data{1};
fclose(fid);
%% Process Demographic responses:
try pp_age=mytab.participant_age(1);
    pp_gender = mytab.participant_gender(1);
catch %stored in responses (later versions)
        disp('ppant age and gender not recovered');
        pp_age=nan;
        pp_gender={'nan'};
end




%% Process PRS-11 Responses: 
% Initialize variables for DASS-10 q's

prs_table = table();
%%
% Initialize output arrays
numRows = length(data);

% Process each data row
Attncounter= 0; % increment for Attn checks. 
Attnperformance= []; % increment for Attn checks. 

for i = 1:numRows
    % Parse this row (if not empty)
    if isempty(data{i})
        disp(['Skipping empty data entry for row ' num2str(i)]);
    else
        rowData = textscan(data{i}, '%s', 'Delimiter', ',');
        rowData = rowData{1};
        

        %% extract the demographics if it didnt already exist.
        if any(contains(rowData,'demographics'))
            
            %extract data: 
            startQ= find(contains(rowData, 'age')); % first entry
            endQ = find(contains(rowData, 'gender')); % first entry
            for iR = startQ:endQ
        
                %extract the number.
                stringR = rowData{iR};
                %replace double quotes with singles temporarily:
                 str_cleaned = strrep(stringR, '""', '''');
                 % enclose:
                 str_cleaned = strrep(str_cleaned, '''', '"');
                 % Now extract all field-value pairs using regular expressions
                 if contains(str_cleaned, 'gender') 
                     % pattern should search for two strings:
                     pattern = '"([^"]+)":"([^"]+)"';                     

                 elseif contains(str_cleaned, 'age'); 
                     pattern = '{?"([^"]+)":"([^"]+)"';
                 end
                 matches = regexp(str_cleaned, pattern, 'tokens');


                 % Process each match
                 for i = 1:length(matches)
                     field_name = matches{i}{1};   % Field name (e.g., "FA_12")
                     
                     if contains(str_cleaned, 'gender') 
                         value = matches{i}{2};  % string value
                         pp_gender = value;
                     elseif contains(str_cleaned, 'age')
                         value = str2double(matches{i}{2}) ; % Numeric value (e.g., 6)
                         pp_age = value;
                     end
                     


                     % store and score
                     if ischar(value)
                        % For string values, create a string array column
                        prs_table.(field_name) = strings(height(prs_table), 1);
                        prs_table.(field_name)= value;
                     else
                        prs_table.(field_name)= value;
                         
                     end

                     end
                 

            end % each response... 
        end % end demographic extraction:

        %% check if we have attn check on this row: 
        if any(contains(rowData,'attncheck'))
        
            Attncounter=Attncounter+1; 

            % was it correct?
            if contains(rowData{attnCheckIdx}, 'true');
                Attnperformance(Attncounter)= 1;
            else
                Attnperformance(Attncounter)= 0;
            end
            
        end

        %% Check if this is the dass10 row
        if  length(rowData)> taskColIdx && any(contains(rowData,'dass10'))
            
            % the Q's are always in order, so find the start and end
            % indices of the PRS questions:
            startQ= find(contains(rowData, ['felt_panic']));
            endQ= find(contains(rowData, 'EOL'));

            for iR = startQ:endQ
        
                %extract the number.
                stringR = rowData{iR};
                %replace double quotes with singles temporarily:
                 str_cleaned = strrep(stringR, '""', '''');
                 % enclose:
                 str_cleaned = strrep(str_cleaned, '''', '"');
                 % Now extract all field-value pairs using regular expressions
                 pattern = '"([^"]+)":(\d+)';
                 matches = regexp(str_cleaned, pattern, 'tokens');

                 % Process each match
                 for i = 1:length(matches)
                     field_name = matches{i}{1};   % Field name (e.g., "FA_12")
                     value = str2double(matches{i}{2});  % Numeric value (e.g., 6)

                     fprintf('Field: %s, Value: %d\n', field_name, value);


                     % store and score
                     prs_table.(['DASS-' field_name])= value;
                 end
                 

            end % each response in the DASS
           
            
        end

        %% here we also have the GAAIS  (some pparticipants).
        if  length(rowData)> taskColIdx && any(contains(rowData,'GAAIS'))
            
            % the Q's are actually randomized! , so find the start and end
            % indices of the PRS questions:
            startQ= find(contains(rowData, '"{"'));
            endQ= find(contains(rowData, '}"'));

            for iR = startQ:endQ
        
                %extract the number.
                stringR = rowData{iR};
                %replace double quotes with singles temporarily:
                 str_cleaned = strrep(stringR, '""', '''');
                 % enclose:
                 str_cleaned = strrep(str_cleaned, '''', '"');
                 % Now extract all field-value pairs using regular expressions
                 pattern = '"([^"]+)":(\d+)';
                 matches = regexp(str_cleaned, pattern, 'tokens');

                 % Process each match
                 for i = 1:length(matches)
                     field_name = matches{i}{1};   % Field name (e.g., "FA_12")
                     value = str2double(matches{i}{2});  % Numeric value (e.g., 6)

                     fprintf('Field: %s, Value: %d\n', field_name, value);


                     % store and score
                     prs_table.(['GAAIS-' field_name])= value;
                 end
                 

            end % each response... 
           
            
        end
    end % empty
end % each row.
% we should now have the DASS-10 and GAAIS results compiled.

%% some basic scoring: 
%https://ebchelp.blueprint.ai/en/articles/9912641-depression-anxiety-and-stress-scale-10-dass-10#h_2f165842c2
% ~TO DO: DASS-10 scored as follows: 
% Q'S basically summed togehter, none reverse scored. 
% 1-felt_panic (A-S)
% 2-no_initiative (D)
% 3-felt_blue (D)
% 4-was_intolerant (A-S)
% 5-nothing_forward (D)
% 6-felt_scared (A-S)
% 7-over_react (A-S)
% 8-felt_worry (A-S)
% 9-no_relax (A-S)
% 10-no_positive (D)

% 11-felt_annoyed
% 12-EOL_thoughts

%NB not all participants have the DASS:
if any(contains(prs_table.Properties.VariableNames, 'DASS'))
    % Depression:
    depression_subscore = nansum([prs_table.("DASS-no_initiative"),...
        prs_table.("DASS-nothing_forward"),...)
        prs_table.("DASS-no_positive")]);

    % Anxiety-Stress:
    anxstress_subscore =  nansum([prs_table.("DASS-felt_panic"),...
        prs_table.("DASS-was_intolerant"),...)
        prs_table.("DASS-felt_scared"),...]);
        prs_table.("DASS-over_react"),...
        prs_table.("DASS-felt_worry"),...
        prs_table.("DASS-no_relax")]);

    % add to table:
    prs_table.DASS_depression_subscore = depression_subscore;
    prs_table.DASS_anxstress_subscore = anxstress_subscore;

else% save as NAN.
prs_table.DASS_depression_subscore = nan;
prs_table.DASS_anxstress_subscore = nan;
end
%% now score the GAAIS: from Shepman & Rodway: https://doi.org/10.1080/10447318.2022.2085400
%20 item score, with positive and negative components. negative are
%negativelt scored.
% 
% Attncheck
% to pass, the response must be "4" (for strongly agree).

if any(contains(prs_table.Properties.VariableNames, 'GAAIS'))
    if prs_table.("GAAIS-A")==4
        GAAIS_attnpassed= 1;
    else
        GAAIS_attnpassed= 0;
    end
    
    %PosAI - scale mean
    posAI_subscore = nanmean([prs_table.("GAAIS-pos_1"),...
        prs_table.("GAAIS-pos_2"),...
        prs_table.("GAAIS-pos_4"),...
        prs_table.("GAAIS-pos_5"),...
        prs_table.("GAAIS-pos_7"),...
        prs_table.("GAAIS-pos_11"),...
        prs_table.("GAAIS-pos_12"),...
        prs_table.("GAAIS-pos_13"),...
        prs_table.("GAAIS-pos_14"),...
        prs_table.("GAAIS-pos_16"),...
        prs_table.("GAAIS-pos_17"),...
        prs_table.("GAAIS-pos_18")]);
    
    %negative AI attitudes - scale mean:
    negAI_subscore = abs(nanmean([prs_table.("GAAIS-neg_3"),...
        prs_table.("GAAIS-neg_6"),...
        prs_table.("GAAIS-neg_8"),...
        prs_table.("GAAIS-neg_9"),...
        prs_table.("GAAIS-neg_10"),...
        prs_table.("GAAIS-neg_15"),...
        prs_table.("GAAIS-neg_19"),...
        prs_table.("GAAIS-neg_20")]));
    
    % add to table
    prs_table.GAAIS_posAI_subscale= posAI_subscore;
    prs_table.GAAIS_negAI_subscale= negAI_subscore;
    prs_table.GAAIS_attnpassed= GAAIS_attnpassed;
else
    % add to table
    prs_table.GAAIS_posAI_subscale= nan;
    prs_table.GAAIS_negAI_subscale= nan;
    prs_table.GAAIS_attnpassed= nan;
end


%% now complete the mytab wrangling:
  %% correct for the bug (no correct_1 or correct_2 - some datasets.).
        allr1 = find(contains(mytab.stimulus, 'first response'));
        allr2 = find(contains(mytab.stimulus, 'second response'));
        allnum = find(contains(mytab.stimulus, 'Dot display')); % numerosity (dot) trials
        allconsol = find(contains(mytab.stimulus, 'consolidated'));
        

        if ~any(contains(mytab.Properties.VariableNames, 'correct_1'))
                        
                disp(['warning! no correct_1 in current table, adding..'])
                % now for each row, add the correct_1 entry:
                
                for itrial = 1:length(allr1) % this is searching all second response.
                    rowindx = allr1(itrial);
                    selectedside = mytab.selected_side{rowindx}; %'left', or 'right'

                    % find previous numerosity trial and what was correct.
                    prevNum = find(contains(mytab.stimulus(1:rowindx), 'Dot display'), 1, "last");
                    % on that trial, which side had more?
                    if mytab.leftDots(prevNum) > mytab.rightDots(prevNum)
                        corside = 'left';
                    elseif mytab.leftDots(prevNum) < mytab.rightDots(prevNum)
                        corside = 'right';
                    else
                        error(['check data on row indx ' num2str(prevNum)]);

                    end

                    if contains(selectedside,corside);
                        thiswas  = 'true';
                    else
                        thiswas = 'false';
                    end

                    mytab.correct_1{rowindx}= thiswas;
                    % also store RT.
                    mytab.rt_1(rowindx) = str2double(mytab.rt{rowindx});

                    %(converted to binary below)

                    %also update the consolidated row entry:
                    cindx = find(allconsol>rowindx,1,'first');
                    crow= allconsol(cindx);

                    %update:

                    mytab.correct_1{crow} =thiswas;
                    mytab.rt_1(crow) =mytab.rt_1(rowindx);

                end



            
            mytab = movevars(mytab, 'correct_1', 'After', 'correct');
            mytab = movevars(mytab, 'rt_1', 'After', 'correct_1');
        end

        % some defensive checks to see if the data in correct_2 exists, or
        % is weird: (unfortunately varies per ppant due to pilotting).

        corr2flag=0;
        if  ~any(contains(mytab.Properties.VariableNames, 'correct_2'))
            corr2flag=1;
        end
        try allnan= isempty(find(~isnan(mytab.correct_2),1));
            if allnan
            corr2flag=1;
            end
        catch
        end

        %
        if corr2flag==1    
        % strip if it exists:
            
            if any(contains(mytab.Properties.VariableNames, 'correct_2'))
                mytab.correct_2=[];
            end
            
            %% continue with this step only if needed (data-md)
            for itrial = 1:length(allr2) % this is searching all second response.
                rowindx = allr2(itrial);
                selectedside = mytab.selected_side{rowindx}; %'left', or 'right'
                
                % find previous numerosity trial and what was correct.
                prevNum = find(contains(mytab.stimulus(1:rowindx), 'Dot display'), 1, "last");
                % on that trial, which side had more?
                if mytab.leftDots(prevNum) > mytab.rightDots(prevNum)
                    corside = 'left';
                elseif mytab.leftDots(prevNum) < mytab.rightDots(prevNum)
                    corside = 'right';
                else
                    error(['check data on row indx ' num2str(prevNum)]);

                end

                if contains(selectedside,corside);
                    thiswas  = 'true';
                else
                    thiswas = 'false';
                end

                mytab.correct_2{rowindx}= thiswas;
                mytab.rt_2(rowindx)= str2double(mytab.rt{rowindx});
                %(converted to binary below)

                %also update the consolidated row entry:
                cindx = find(allconsol>rowindx,1,'first');
                crow= allconsol(cindx);

                %update:

                mytab.correct_2{crow} =thiswas;
                mytab.rt_2(crow) =mytab.rt_2(rowindx);

            end
        end

%%
            mytab = movevars(mytab, 'correct_2', 'After', 'selectedSide_second');
            mytab = movevars(mytab, 'rt_2', 'After', 'correct_2');
        %% convert to binary most true-false columns:
        forcols ={'correct', 'correct_1', 'adviceCorrect', 'correct_2'};
        for icol=1:length(forcols);

            % extract the data.
            tmpdata = nan(height(mytab),1);
            for itrial = 1:height(mytab);
                tmp= mytab.([forcols{icol}])(itrial);
                if isempty(tmp{1})
                    tmpdata(itrial,1)= nan;
                elseif contains(tmp{1}, 'true');
                    tmpdata(itrial,1)= 1;
                elseif contains(tmp{1}, 'false');
                    tmpdata(itrial,1)= 0;

                end
            end
            % rename and save after relevant data.:
            mytab.([forcols{icol} '_binary'])= tmpdata;
            mytab= movevars(mytab, [forcols{icol} '_binary'], 'after', forcols{icol});

        end


        %% adjust confidence values before zscoring.

        % we are mainly interested in the consolidated data rows.
        consolrows = find(contains(mytab.task, 'consolidated'));
        for irow = 1: length(consolrows);

            %for each row, subtract 50 from confidence values:
            userow = consolrows(irow);
            mytab.confidence_first(userow) = mytab.confidence_first(userow)-50;
            mytab.confidence_second(userow) = mytab.confidence_second(userow)-50;
            % adjust to be - for left, + for right?
            % if contains(mytab.selectedSide_first(userow),'left');
            %     mytab.confidence_first(userow) = mytab.confidence_first(userow)*-1;
            % end
            % if contains(mytab.selectedSide_second(userow),'left');
            %     mytab.confidence_second(userow) = mytab.confidence_second(userow)*-1;
            % end

        end

        %% now zscore all.
        allc1 = mytab.confidence_first;
        allc2 = mytab.confidence_second;
        allconf= [allc1;allc2];
        % not zscore doesnt work with nan.

        confreal = allconf(~isnan(allconf));
        zconf = zscore(confreal);
        %%
        % now the tricky bit, for each original conf value, find the corresponding
        % zscore and place in the table.
        % only use the consol rows.
        % first create nan rows:
        mytab.confidence_first_z=nan(height(mytab),1);
        mytab.confidence_second_z=nan(height(mytab),1);
        %

        for itrial=1:length(consolrows)
            ir = consolrows(itrial);
            ogC1= mytab.confidence_first(ir);
            ogC2= mytab.confidence_second(ir);

            % find the index in confreal.
            indxs = dsearchn(confreal, [ogC1, ogC2]');
            %use these indexes for the zscores:
            znow = zconf(indxs,1);
            mytab.confidence_first_z(ir) = znow(1);
            mytab.confidence_second_z(ir) = znow(2);
        end
        % when all said and done, move the zscored cols to the right spot.
        mytab = movevars(mytab, 'confidence_first_z', 'After', 'confidence_first');
        mytab = movevars(mytab, 'confidence_second_z', 'After', 'confidence_second');

        % include change in conf.
        confChange = mytab.confidence_second_z - mytab.confidence_first_z;
        mytab.confChange_z = confChange;
        mytab = movevars(mytab, 'confChange_z', 'After', 'confidence_second_z');


        % also include whether this change in confidence is consistent with
        % decision?? 

        %%
      

%% SAVE 
        disp('saving grand participant table')...
        % with mytab now clean, lets save for reuse: 
        save(['Participant_' pnum], 'mytab', ...
            'prs_table', ...
            'Attncounter', ...
            'Attnperformance', 'pnum');

        
        
        %% now for condensed data from the experiment.
        disp('creating condensed participant data table...')

        consolidated_rows = find(contains(mytab.stimulus,'consolidated-trial-data'));
        ctable = mytab(consolidated_rows,:);

        %%
        % need to do some reordering.
        % delete the bloat columns, single entries we dont need.
        del_columns={'trial_type',...
            'trial_index',...
            'participant_age',...
            'participant_gender',...
            'correct_binary',... % not in ctable
            'attention_check_correct',...
            'target_button', ...
            'question_order',...
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
        % choice selected from the previous '1' (and advice direction).

        % change adviceChoice to nan for adviceConds 2,3.
        % change adviceCorrect to nan for adviceConds 2,3.
        % change adviceDir to nan for adviceCond 3 (no advice)

        for itrial= 1:height(ctable)
            
            % if they didnt choose to see advice, adjust table:
            if ctable.adviceCond(itrial)==1 && contains(ctable.adviceChoice{itrial}, 'no');
                % if they didnt choose to see advice, adjust table:
                ctable.adviceDir{itrial}= 'none';
                ctable.adviceCongruent(itrial)= nan;
                ctable.adviceCorrect{itrial}= 'nan';

            end
            % if it was a forced choice trial (2), or no advice trial (3)
            if ctable.adviceCond(itrial)>1
                ctable.adviceChoice{itrial}='NaN'; % text to make string search easier
                % extra step:
                if ctable.adviceCond(itrial)==3   %(no advice)  
                    ctable.adviceDir{itrial}= 'none';
                    ctable.adviceCorrect{itrial}= 'nan';
                    ctable.adviceCongruent(itrial)= nan;
                end
            end
        end
        %% now correct the binary columns based on advice:
        
        nanrows = find(contains(ctable.adviceCorrect, 'nan'));
        ctable.adviceCorrect_binary(nanrows)=nan;
    %% include a column for whether advice was congruent with decision 1.
        adviceCongr=[];
        for itrial= 1:height(ctable)
            dir1 = ctable.selectedSide_first(itrial);

            if contains(ctable.adviceDir(itrial), 'none');
                adviceCongr(itrial,1)=nan;
            
            elseif contains(ctable.adviceDir(itrial),dir1);
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

%% now we can reclassify the change in confidence, based on whether 
% the direction of the original decision. 
% e.g. (+) for increased conf in original decision
% e.g. (-) for decreased conf in original decision

%pseudo:
%- find all advice offered trials
%- for 'left' on decision 1, multiply by -1

%
allconfchange_z = ctable.confChange_z;
%when was the first decision on left?
left1decision = find(contains(ctable.selectedSide_first, 'left')); 
% for all these ros, a negative change in confidence is an increase in
% decision. and vice versa. 
allconfchange_z(left1decision) = allconfchange_z(left1decision).*-1 ;
%inlcude in table. 
ctable.confChangevsFirst_z=  allconfchange_z;
ctable= movevars(ctable, 'confChangevsFirst_z', 'After', 'confChange_z');


%% similarly, include a change of mind column.
firstSide = ctable.selectedSide_first;
secondSide = ctable.selectedSide_second;
COM = strcmp(firstSide,secondSide);
ctable.changeOfMind = abs(COM-1); % since 1 was a match.
ctable= movevars(ctable, 'changeOfMind', 'After', 'confChangevsFirst_z');
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

        %% SAVE now save this revised condensed table.
        disp('Saving condensed table');
        % with mytab now clean, lets save for reuse: 
        save(['Participant_' pnum], 'ctable', '-append');

 %% also prep a participant data structure for easy ouput later.
    % things we want to save include: 
    % demographics, survey results, attn check performance. etc
    % propY selected (total)

    % Accuracy on C1 (total)
    % Accuracy on C2 (total)
    
    % then various DVs we can split by subgrops:
    % Acc on C1
    % RT on C1 
    % Conf on C1
    
    % Acc on C2
    % RT on C2 
    % Conf on C2 
       
    % Change in Conf (relative to first decision)    
    
    ParticipantData=[];
    
    %% First, extract demongraphcis and survey responses:
    ParticipantData.pp_Age= pp_age;
    ParticipantData.pp_gender= {pp_gender};

    % add DASS subscores here.
    ParticipantData.DASS_depression= prs_table.DASS_depression_subscore;
    ParticipantData.DASS_anx= prs_table.DASS_anxstress_subscore;

    %GAAIS:
    ParticipantData.GAAIS_attnPassed= prs_table.GAAIS_attnpassed;
    ParticipantData.GAAIS_posAI= prs_table.GAAIS_posAI_subscale;  
    ParticipantData.GAAIS_negAI= prs_table.GAAIS_negAI_subscale;

    % how many attn checks in exp?
    ParticipantData.pp_AttnPerformance= mean(Attnperformance);
    ParticipantData.pp_nAttnChecks= length(Attnperformance);
    
    % also calculate overall likelihood of seeing advice:
    %prop selected
        choiceTrials= find(ctable.adviceCond==1);
        choices = ctable.adviceChoice(choiceTrials);
        choseY = find(contains(choices, 'yes'));        
        propY= length(choseY)/length(choiceTrials);
    
    ParticipantData.proportionYesSelected = propY;

%%
    %overall performance?
    ParticipantData.choice1_accuracy = nanmean(ctable.correct_1_binary);
    ParticipantData.choice2_accuracy = nanmean(ctable.correct_2_binary);
    
    
    %next we will split various DVs by subgroups:
    
    %create trial lists for data aggregation:
    
    correct_first= find(ctable.correct_1_binary==1);
    error_first= find(ctable.correct_1_binary==0);
    
    correct_second = find(ctable.correct_2_binary==1);
    error_second = find(ctable.correct_2_binary==0);
        
    choiceY_idx = intersect(find(ctable.adviceCond==1), find(contains(ctable.adviceChoice,'yes')));
    choiceN_idx = intersect(find(ctable.adviceCond==1), find(contains(ctable.adviceChoice,'no')));
    forced_idx = find(ctable.adviceCond==2);
    noAdv_idx = find(ctable.adviceCond==3);

    %check indexing:
    advCongruent = find(ctable.adviceCongruent==1);
    advIncongruent = find(ctable.adviceCongruent==0);

    % advice x congruence with c1 
    choiceY_congr = intersect(choiceY_idx, advCongruent );
    choiceY_incongr = intersect(choiceY_idx, advIncongruent );
    
    forced_congr = intersect(forced_idx, advCongruent);
    forced_incongr = intersect(forced_idx, advIncongruent); 


    searchList = {correct_first,error_first,...
        correct_second, error_second,...
        choiceY_idx, choiceN_idx, forced_idx, noAdv_idx,...
        advCongruent, advIncongruent, ...
        choiceY_congr, choiceY_incongr,...
        forced_congr, forced_incongr};

    %%
    % create descriptors for the field names/ columns. 
    flist ={'correct_first', 'error_first',...
        'correct_second','error_second',...
        'choseYesAdvice','choseNoAdvice','forcedYesAdvice','NoAdvice',...}
        'allAdviceCongruent','allAdviceIncongruent',...}
        'choseYesAdvice_Congr', 'choseYesAdvice_Incongr',...
        'forcedYesAdvice_Congr', 'forcedYesAdvice_Incongr'};
        


    for ilist = 1:length(searchList)
        currlist = searchList{ilist};
        currLabel = flist{ilist};

%First store the data for choice 1: 
%Acc,RT,conf on first choice:
ParticipantData.(['choice1_accuracy_' currLabel ]) = nanmean(ctable.correct_1_binary(currlist));
ParticipantData.(['choice1_rt_' currLabel ]) = nanmean(ctable.rt_1(currlist));
ParticipantData.(['choice1_abs_conf_' currLabel ]) = nanmean(abs(ctable.confidence_first_z(currlist)));

%Acc,RT,conf on second choice:
ParticipantData.(['choice2_accuracy_' currLabel ]) = nanmean(ctable.correct_2_binary(currlist));
ParticipantData.(['choice2_rt_' currLabel ]) = nanmean(ctable.rt_2(currlist));
ParticipantData.(['choice2_abs_conf_' currLabel ]) = nanmean(abs(ctable.confidence_second_z(currlist)));

%changes in CONF:
ParticipantData.(['confChange_abs_z_' currLabel]) = nanmean(ctable.confChangevsFirst_z(currlist)); % dont take absolute, because here the sign has meaning
% COM likelihood:
ntotal = length(currlist);
nhappened = nansum(ctable.changeOfMind(currlist));
ParticipantData.(['changeOfMind_prop_' currLabel]) = nhappened/ntotal;

    end

    % include demographics, attention checks, and survey responses: 
    %%
% save processed data:
    save(['Participant_' pnum], 'ParticipantData', '-append');

    
    end  % per ppant
end % job import and wrangle


%% Participant level plots:


if jobs.plot_PP_summary
    
%% summary can be a single figure, 
cd(datadir)
ppantfols = dir([pwd filesep 'Participant_*']);
%

for ippant=1%:length(pfols);
cd(datadir)
load(ppantfols(ippant).name);

set(gcf,'color','w', 'units','normalized','position', [0 0 1 1]); clf
shg
%

%% first subplot, overall performance on C1 and C2.
totalAcc= [ParticipantData.choice1_accuracy, ParticipantData.choice2_accuracy];
subplot(3,3,1);
bar(1:2, totalAcc);
title('Total Accuracy')
ylabel('Proportion correct');
set(gca,'XTickLabels', {'first choice', 'second choice'}, 'fontsize', 15);

%% next the RTs when Correct and Error, for first
dataRT = [ParticipantData.choice1_rt_correct_first,ParticipantData.choice1_rt_error_first];
subplot(3,3,2);
bar(1:2, dataRT);
title('First choice RT')
% legend({'correct', 'error'})
ylabel('ms');
set(gca,'XTickLabels', {'correct', 'error'}, 'fontsize', 15);

%% conf on first choice

dataConf = [ParticipantData.choice1_abs_conf_correct_first,ParticipantData.choice1_abs_conf_error_first];
subplot(3,3,3);
bar(1:2, dataConf);
title('First choice Conf ')
% legend({'correct', 'error'})
ylabel('zscored Conf');
set(gca,'XTickLabels', {'correct', 'error'}, 'fontsize', 15);
%% Final Accuracy by advice:
DVsare = {'choice2_accuracy_', 'confChange_abs_z_', 'changeOfMind_prop_'};
titlesare= {'Final Accuracy', 'Change in Confidence', 'Changes of mind'};
ylabelsare= {'Proportion correct', 'zscored Conf', 'COM prop.'};

adviceLabels= {'choseYesAdvice','choseNoAdvice', 'forcedYesAdvice','NoAdvice'};
for iDV= 1:length(DVsare)

dataAdv = [ParticipantData.([DVsare{iDV} 'choseYesAdvice']),...
    ParticipantData.([DVsare{iDV} 'choseNoAdvice']),...
    ParticipantData.([DVsare{iDV} 'forcedYesAdvice']),...
    ParticipantData.([DVsare{iDV}  'NoAdvice'])];


subplot(3,3,iDV+3);
bar(1:length(dataAdv), dataAdv);
title(titlesare{iDV})
ylabel(ylabelsare{iDV});
set(gca,'fontsize', 15, 'XtickLabel', adviceLabels);
end

%% now split by congruent or not?
adviceLabels= {'choseYes','forcedYes'};
for iDV= 1:length(DVsare)

dataAdv = [ParticipantData.([DVsare{iDV} 'choseYesAdvice_Congr']),...
    ParticipantData.([DVsare{iDV} 'forcedYesAdvice_Congr']);...
    ParticipantData.([DVsare{iDV} 'choseYesAdvice_Incongr']),...
    ParticipantData.([DVsare{iDV}  'forcedYesAdvice_Incongr'])];


subplot(3,3,iDV+6);
bar(1:2, dataAdv);
title(titlesare{iDV})
ylabel(ylabelsare{iDV});
legend('advice Congruent', 'advice Incogruent')
set(gca,'fontsize', 15, 'XtickLabel', adviceLabels);
end


%%
cd(datadir)
cd ../Figures

print('-dpng', ['Participant_' pnum '_summary'])

end

end

if jobs.createGFX_table
   %% 
   cd(datadir)
    

    ppantfols = dir([pwd filesep 'Participant_*']);
%%
    GFX_table=[];

for ippant=1:length(pfols);
cd(datadir)
load(ppantfols(ippant).name);
%conver tto table and append. 
ptable = struct2table(ParticipantData);

% remove the attn checls
if ippant==1
    GFX_table = ptable;
else
    GFX_table= [GFX_table; ptable];
end

    
end % end ppant

% here also add some ppants (random) based on jitter.
% thinking mean +- SD

GFX_table_Pseudo = GFX_table;

nCreate = 10; 

for idata= 1:nCreate
% take the average of real (per field), and then add some jitter.
    
    pnew= table(); % new per ppant.
    vrnames = ptable.Properties.VariableNames;
    for ifield = 1:length(vrnames);
        if ifield~=2 % avoid the gender column.
        tmpd= nanmean(GFX_table.([vrnames{ifield}]));
        %lowerbound: a
    a = tmpd - 2*nanstd(GFX_table.([vrnames{ifield}]));
    % upper bound: b
    b = tmpd + 2*nanstd(GFX_table.([vrnames{ifield}]));

    %generate singly number between range (a - b)
    % Generate a single random number
    random_value = a + (b-a) * rand();
        else
            random_value = nan;
        end
    % add this jittered value (somewhere within -2SD + 2SD of mean), to
    % table:
    pnew.([vrnames{ifield}]) = random_value;
    
    end

    % now add to GFXtable:
    GFX_table_Pseudo = [GFX_table_Pseudo; pnew];

end

%% export
writetable(GFX_table, 'GFX_table_Pseudo.csv')
writetable(GFX_table, 'GFX_table.csv')
end

%% now that we have the GFX_table, include group level plots.
if jobs.plotGFX==1
%% following the format above, nowusing GFX_table:

set(gcf,'color','w', 'units','normalized','position', [0 0 1 1]); clf
shg
%

%% first subplot, overall performance on C1 and C2.
totalAcc= [GFX_table.choice1_accuracy, GFX_table.choice2_accuracy];
subplot(3,3,1);
bar(1:2, mean(totalAcc,1));
hold on;
errorbar(1:2,mean(totalAcc,1), CousineauSEM(totalAcc), 'color','k','LineWidth',2, 'LineStyle','none')

title('Total Accuracy')
ylabel('Proportion correct');
set(gca,'XTickLabels', {'first choice', 'second choice'}, 'fontsize', 15);

%% next the RTs when Correct and Error, for first and second (combined?)
dataRT = [GFX_table.choice1_rt_correct_first,GFX_table.choice1_rt_error_first];


subplot(3,3,2);
bar(1:2, nanmean(dataRT,1));
hold on;
errorbar(1:2,nanmean(dataRT,1), CousineauSEM(dataRT), 'color','k','LineWidth',2, 'LineStyle','none')
title('First choice RT')
% legend({'correct', 'error'})
ylabel('ms');
set(gca,'XTickLabels', {'correct', 'error'}, 'fontsize', 15);

%% conf on first choice

dataConf = [GFX_table.choice1_abs_conf_correct_first,GFX_table.choice1_abs_conf_error_first];
subplot(3,3,3); cla
bar(1:2, mean(dataConf)); 
hold on;

errorbar(1:2, mean(dataConf,1), CousineauSEM(dataConf), 'color','k','LineWidth',2, 'LineStyle','none')

title('First choice Conf ')
ylabel('zscored Conf');
set(gca,'XTickLabels', {'correct', 'error'}, 'fontsize', 15);
%% Final Accuracy by advice:
DVsare = {'choice2_accuracy_', 'confChange_abs_z_', 'changeOfMind_prop_'};
titlesare= {'Final Accuracy', 'Change in Confidence', 'Changes of mind'};
ylabelsare= {'Proportion correct', 'zscored Conf', 'COM prop.'};

adviceLabels= {'choseYesAdvice','choseNoAdvice', 'forcedYesAdvice','NoAdvice'};
for iDV= 1:length(DVsare)

dataAdv = [GFX_table.([DVsare{iDV} 'choseYesAdvice']),...
    GFX_table.([DVsare{iDV} 'choseNoAdvice']),...
    GFX_table.([DVsare{iDV} 'forcedYesAdvice']),...
    GFX_table.([DVsare{iDV}  'NoAdvice'])];


subplot(3,3,iDV+3);
bar(1:size(dataAdv,2), nanmean(dataAdv,1));
hold on;
errorbar(1:size(dataAdv,2),nanmean(dataAdv,1), CousineauSEM(dataAdv), 'color','k','LineWidth',2, 'LineStyle','none')
% errorbar_groupedfit(mean(dataAdv,1), CousineauSEM(dataAdv));
title(titlesare{iDV})

ylabel(ylabelsare{iDV});
set(gca,'fontsize', 15, 'XtickLabel', adviceLabels);
end

%% now split by congruent or not?
adviceLabels= {'choseYes','forcedYes'};
for iDV= 1:length(DVsare)

dataAdv = [GFX_table.([DVsare{iDV} 'choseYesAdvice_Congr']),...
    GFX_table.([DVsare{iDV} 'forcedYesAdvice_Congr']),...
    GFX_table.([DVsare{iDV} 'choseYesAdvice_Incongr']),...
    GFX_table.([DVsare{iDV}  'forcedYesAdvice_Incongr'])];

mdata = nanmean(dataAdv,1);
semdata = CousineauSEM(dataAdv);
plotmdata= [mdata(1), mdata(2); mdata(2), mdata(4)];
plotsemdata= [semdata(1), semdata(2); semdata(2), semdata(4)];
subplot(3,3,iDV+6);
bar(1:2, plotmdata);
hold on;
errorbar_groupedfit(plotmdata, plotsemdata);
% errorbar(1:2,mean(dataAdv,1), CousineauSEM(dataAdv), 'color','k','LineWidth',2, 'LineStyle','none')
title(titlesare{iDV})
ylabel(ylabelsare{iDV});
legend('advice Congruent', 'advice Incogruent')
set(gca,'fontsize', 15, 'XtickLabel', adviceLabels);
end


%%
cd(datadir)
cd ../Figures

print('-dpng', ['GFX_n=' num2str(nppants) '_summary'])


end

% 
% %%%%%%%
% if jobs.plot_ppPractice ==1
% 
%     %% first sanity check, staircase performance during practice
% 
%         exp_startrow = find(contains(mytab.stimulus,'expInstruc'), 1,'last');
%         % only one trial type in practice (confidence_first)
%         conf1_rows= find(contains(mytab.task(1:exp_startrow), 'confidence_first'));
%         %numerosity rows have the dot diff and staircase performance.
%         numer_rows = find(contains(mytab.task(1:exp_startrow) , 'numerosity'));
%         numertab = mytab(numer_rows,:);
% 
%         % now plot for sanity check:
%         clf;
%         ntrials = 1:height(numertab);
%         subplot(221);
%         plot(ntrials, numertab.difference)
%         ylabel('dotdifference per trial');
%         title('Practice')
%         shg
% 
%         % rolling average
%         estAv = nan(length(ntrials),1);
%         conftab = mytab(conf1_rows,:);
%         for itrial = ntrials
%             estAv(itrial,1)= sum(conftab.correct_binary(1:itrial))/ itrial;
%         end
%         %
%         subplot(222);
%         plot(ntrials, estAv, 'r-o')
%         ylabel('rolling accuracy');
%         title('Practice')
%         shg
% end
% %%
% if jobs.plot_ppStaircase==1
%      %% plot dot diff and staircase:
% 
%         clf;
%         ntrials = 1:height(ctable);
%         subplot(221);
%         plot(ntrials, ctable.dotsDiff);
%         ylabel('dotdifference per trial');
%         title('Experiment')
%         shg
%         % each bloc had 15 trials, so mark block start end. 
%         hold on;
%         imark = 15:15:225;
%         for im=  1:length(imark)
%         plot([imark(im), imark(im)], [ylim], 'r:')
% 
%         end
% 
%         % rolling average
%         [estAv, estAv2] = deal(nan(length(ntrials),1));
% 
%         for itrial = ntrials
%             estAv(itrial,1)= nansum(ctable.correct_1binary(1:itrial))/ itrial;
%             estAv2(itrial,1)= nansum(ctable.correct_2binary(1:itrial))/ itrial;
%         end
% 
%         %
%         subplot(222);
%         plot(ntrials, estAv, 'b-o'); hold on;
%         plot(ntrials, estAv2, 'r-o'); hold on;
%         ylabel('rolling accuracy');
%         title('Experiment')
%         legend('first resp', 'second');
%         shg
% end
% 
% if jobs.plot_ppAccandConf_basic
%     %% what is the total accuracy on c1, then c2 by option?
%         accC1 = sum(ctable.correct_1binary)/ height(ctable);
%         accC2 = sum(ctable.correct_2binary)/ height(ctable);
% 
% 
%         choiceY_idx = intersect(find(ctable.adviceCond==1), find(contains(ctable.adviceChoice,'yes')));
%         choiceN_idx = intersect(find(ctable.adviceCond==1), find(contains(ctable.adviceChoice,'no')));
%         forced_idx = find(ctable.adviceCond==2);
%         noAdv_idx = find(ctable.adviceCond==3);
% 
%         subplot(2,2,3);
%         bar(1:2, [accC1, accC2]);
%         ylabel('Accuract on first and second choice');
%         ylabel('proportion correct');
%         set(gca,'xticklabels', {'first', 'second'});
%         xlabel('response order');
%         ylim([.5 1]);
% 
%         %acc per type
%         chY_acc = sum(ctable.correct_2binary(choiceY_idx))/ length(choiceY_idx);
%         chN_acc = sum(ctable.correct_2binary(choiceN_idx))/ length(choiceN_idx);
%         forced_acc= sum(ctable.correct_2binary(forced_idx))/ length(forced_idx);
%         noAdv_acc= sum(ctable.correct_2binary(noAdv_idx))/ length(noAdv_idx);
% 
%         subplot(224);
%         bar(1:4, [chY_acc, chN_acc, forced_acc, noAdv_acc]);
%         xlabel('advice condition');
%         ylabel('proportion correct');
%         set(gca,'xticklabels', {'chose Yes', 'chose No', 'forced', 'no advice'})
%         shg
% 
% 
%         %% does confidence track accuracy?
%         respwas ={'first', 'second'};
%         PFX_Conf=[];
%         PFX_RT=[];
%         for iresp = 1:2
%             corr_idx = find(ctable.(['correct_' num2str(iresp) 'binary'])==1);
%             err_idx = find(ctable.(['correct_' num2str(iresp) 'binary'])==0);
% 
%             % collect responses. - absolute for confidence data.
%             confData = abs(ctable.(['confidence_' respwas{iresp}]));
% 
%             rtData= ctable.(['rt_' num2str(iresp)]);
% 
%             meanConf = [nanmean(confData(corr_idx)),...
%                 nanmean(confData(err_idx))];
% 
%             meanRT = [nanmean(rtData(corr_idx)), ...
%                 nanmean(rtData(err_idx))];
% 
%             PFX_Conf(iresp,1:2) = meanConf;
%             PFX_RT(iresp,1:2) = meanRT;
%         end
%         %%plot
% 
%         subplot(2,2,3);
%         bar(PFX_Conf); ylabel('mean confidence')
%         set(gca,'xticklabels', {'first response', 'second response'});
%         legend('correct', 'error', 'Location', 'northoutside')
% 
%         subplot(2,2,4);
%         bar(PFX_RT);ylabel('mean RTs');
%         % ylim([0 1]);
%         set(gca,'xticklabels', {'first response', 'second response'});
%         legend('correct', 'error', 'Location', 'northoutside')
%         shg
% end
% 
% if jobs.plot_ppAgreegmentwithAdvice    
% %% next figure. plot agreement:
%         % pseudo:
%         % Proportion selected side = advice side (?)
%         % Confidence change (congr vs incongr)
% 
%         %prop selected
%         choiceTrials= find(ctable.adviceCond==1);
%         choices = ctable.adviceChoice(choiceTrials);
%         choseY = find(contains(choices, 'yes'));
%         choseN = find(contains(choices, 'no'));
%         propY= length(choseY)/length(choiceTrials);
%         propN= length(choseN)/length(choiceTrials);
%         % initial Confidence for each.
%         conf1onchoice = abs(ctable.confidence_first(choiceTrials))
%         conf_chY = conf1onchoice(choseY);
%         conf_chN = conf1onchoice(choseN);
%         clf;
%         subplot(1,3,1)
%         bar(1, [propY;propN], 'stacked');
%         legend('yes', 'no');
%         ylabel('proportion'); title('would you like advice?')
% 
%         subplot(132);
%         bar(1:2, [mean(conf_chY), mean(conf_chN)]);
%         % ylim([50 100]);
%         title('initial confidence'); set(gca,'xticklabels', {'yes', 'no'}); xlabel('Would you like advice?');
%         ylabel('initial confidence')
% 
%         % now final confidence for each.
%         conf2onchoice = abs(ctable.confidence_second(choiceTrials));
%         conf_chY = conf2onchoice(choseY);
%         conf_chN = conf2onchoice(choseN);
% 
%         subplot(133);
%         bar(1:2, [mean(conf_chY), mean(conf_chN)]);
%         % ylim([50 100]);
%         title('final confidence'); set(gca,'xticklabels', {'yes', 'no'}); xlabel('Would you like advice?');
%         ylabel('final confidence')
%         shg
% end
% %% next figure. plot agreement:
% % pseudo:
% % Proportion selected side = advice side (?)
% % Confidence change (congr vs incongr)
% 
% %prop selected
% choiceTrials= find(ctable.adviceCond==1);
% choices = ctable.adviceChoice(choiceTrials);
% choseY = find(contains(choices, 'yes'));
% choseN = find(contains(choices, 'no'));
% propY= length(choseY)/length(choiceTrials);
% propN= length(choseN)/length(choiceTrials);
% % initial Confidence for each.
% conf1onchoice = abs(ctable.confidence_first(choiceTrials));
% conf_chY = conf1onchoice(choseY);
% conf_chN = conf1onchoice(choseN);
% 
% 
% clf;
% subplot(1,3,1)
% bar(1, [propY;propN], 'stacked');
% legend('yes', 'no');
% ylabel('proportion'); title('would you like advice?')
% 
% subplot(132);
% bar(1:2, [mean(conf_chY), mean(conf_chN)]);
% % ylim([50 100]);
% title('initial confidence'); set(gca,'xticklabels', {'yes', 'no'}); xlabel('Would you like advice?');
% ylabel('initial confidence')
% 
% % now final confidence for each.
% conf2onchoice = abs(ctable.confidence_second(choiceTrials));
% conf_chY = conf2onchoice(choseY);
% conf_chN = conf2onchoice(choseN);
% 
% subplot(133);
% bar(1:2, [mean(conf_chY), mean(conf_chN)]);
% % ylim([50 100]);
% title('final confidence'); set(gca,'xticklabels', {'yes', 'no'}); xlabel('Would you like advice?');
% ylabel('final confidence')
% shg
% 
% 
