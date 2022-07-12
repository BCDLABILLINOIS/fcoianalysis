%Matlab analytic code that forms the basis for defining fcois, extracting fixed and fcoi data for analysis, running central inferential
%statistics, and calculating consistency in fcoi definitions for simulated
%dataset (i.e., study 3); The same basic code is used for all studies in
%Liu,Sánchez Hernández, Ting, & Hyde: Comparing fixed-array and functionally-defined channel of interest approaches to infant functional near-infrared spectroscopy data

%NOTE: Data plotting and other statisitic/metric extraction for reporting handled
%in other scripts--available upon request from authors(dchyde@illinois.edu)

%% 
close all;
clear all;
%INPUT DESIRED ANALYSIS SETTINGS HERE
topchan=1; %how many channels to include in fcoi definition (=1,2,or 3)
restricttochans=0; %anatomical restriction 1==yes/0==no
clean=1; %if clear out runs with no useable data;=0 if no cleaning/may not work


%%
%%load pre-processed results nirs data Homer2
load('S3groupResults.mat'); 
%%
%other dataset-specific settings--will change with dataset
timest=2; %start of window of interest for analysis (in seconds)
timeend=14;%end of window of interest for analysis (in seconds)
samplingrate=50; %sampling rate in Hz
baseline=2; %baseline before stim onset to include (in seconds)
permis =1; %variable placeholder for future/possible permutation methods---not relevant to current analyses here
%calcuate some additional setttings based on choices above
activeselect = 1;%if you want to look at response to "active condition--1" *would need code change for reverse effect)
topchanrange=1:topchan; %1:X vector with # of channels used for fcoi def.
if restricttochans==0
    whatrestrict = 0;%not used if restricttochans = 0;
elseif restricttochans==1
    whatrestrict = 1:5; %restrict to channels 1:5 if restricttochans = 1;will need to change if probe/dataset changes
end
%futher define time window of interest
t = [timest timeend]*samplingrate; %input window of interest x sampling rate
base = baseline*samplingrate;
tt = t+base;
%find number of runs of data
runz = numel(group.subjs(1,:));
%find subject numbers and number of subjects from names in Homer group file
for i = 1:runz
    sublist(i) = sscanf(group.subjs(i).name,'synfNIRS2chanvar%d');%text change to match subject names in Homer groupresults file
end
subnumbs = (unique(sublist));%find unique subject numbers
numbsubs = numel(subnumbs);%find number of unique subjects

%%
%%some dataset-specific cleaning procedures:
%%pulls mean from window of interest for each line of data (i.e., run) for all channels for all conditions
for eachline = 1:numel(sublist) %subs
        subrundata(eachline,:,1) = squeeze(nanmean(group.subjs(eachline).procResult.dcAvg(tt(1):tt(2),1,:,1)));%active condition
        subrundata(eachline,:,2) = squeeze(nanmean(group.subjs(eachline).procResult.dcAvg(tt(1):tt(2),1,:,2)));%control condition
end
%%%%find all runs where both conditions are NaN/not good
for iv = 1:numel(sublist)%for each run of data
    if isnan(nanmean(subrundata(iv,:,1)))==1 && isnan(nanmean(subrundata(iv,:,2)))==1
        subrundatacleanout(iv,1) = 1;%if all channels/all conditions in run are bad, set up to remove.
    else subrundatacleanout(iv,1) = 0; %otherwise set up to keep
    end
end
%clear runs without any data (NaNs--bad for all conditions and all channels)
if clean ==1
    L=logical(subrundatacleanout);
    subrundatacleaned = subrundata;%set up temp file of subrundata
    sublistcleaned = sublist;%set up temp file of sublist
    subrundatacleaned(L,:,:)=[];%take out bad runs
    sublistcleaned(L)=[];%take out subject indexes corresponding to bad runs
    subrundata=[];%make sure original variable is cleared
    sublist=[];%make sure original variable is cleared
    subrundata=subrundatacleaned; %rename temp variable with original name: subrundata
    sublist =sublistcleaned;%rename temp variable with original name: sublist
end

%%
%anatomical restriction settings
channunber = numel(subrundata(1,:,1));%figure out number of channels in current data
%reduce subrundata if anatomical restriction is in place
if restricttochans==1
    subrundata=subrundata(:,whatrestrict,:);
    channunber = numel(subrundata(1,:,1));
else
end

%%
%%configure data for fcoi analyses
%%pull means for each run for each channel for each subject and then calculates differences and sorts
%%to find most responsive channels
for eachsub = 1:numbsubs
    subchandata(eachsub,:,:) = squeeze(nanmean(subrundata(sublist==subnumbs(eachsub),:,:),1));%averages all runs associated with given subject
    if activeselect==1
        subdiff = subchandata(eachsub,:,1)-subchandata(eachsub,:,2);%calculate difference if examining active>control for all channels
    elseif activeselect==2
        subdiff = subchandata(eachsub,:,2)-subchandata(eachsub,:,1);%calcuate difference if examining control>active for all channels
    end
    [MS, IZ] = sort(subdiff,'descend','MissingPlacement','last'); %sort differences from largest with associated channel numbers
    %MS==difference values
    %IZ = associated channel numbers
    if isnan(MS(1))==1%if top difference is NaN there is a problem with data
        subchanmax(eachsub,:)=NaN;%mark as NaN
    elseif ~isnan(MS(1))==1%if top differences is actual value.
        subchanmax(eachsub,:) = IZ; %record sorted list of channels for given subject
    end
end

%%
%%extract data for fcoi methods
%%between-subjects leave-one-subject-out method
for permz = 1:permis % placeholder for possible/future permutation methods--not relevant to this/==1 here.
    for eachsubjz = 1:numbsubs
        leaveonsubout = subchandata;%make new temp. variable equal to all data %%work with subchandata
        fulldatasetloso = subchandata;
        leaveonsubout(eachsubjz,:,:) = [];%delete target subject from temp variable
        meanleftoutsub = squeeze(nanmean(leaveonsubout,1)); %find means for remaining data for each channel
        if activeselect==1 %if active condition
            meanleftoutsubdiff = meanleftoutsub(:,1)-meanleftoutsub(:,2);%calculate diff.
        elseif activeselect==2 %if control condition
            meanleftoutsubdiff = meanleftoutsub(:,2)-meanleftoutsub(:,1);%calculate diff.
        end
        [ME, IE] = sort(meanleftoutsubdiff,'descend','MissingPlacement','last');
        %pull means for each channel up to number of top channels taken
        for pulltopchan = 1:topchan
            if isnan(ME(pulltopchan,1))==1  %if nan, mark as such
                leftoutsubpullch(permz,eachsubjz,pulltopchan,1)= NaN;
                leftoutsubpullch(permz,eachsubjz,pulltopchan,2)= NaN;
                sCOI(permz,eachsubzz,pulltopchan) = NaN;
            elseif ~isnan(ME(pulltopchan,1))==1 %if not nan, pull data
                leftoutsubpullch(permz,eachsubjz,pulltopchan,1) = fulldatasetloso(eachsubjz,IE(pulltopchan),1);
                leftoutsubpullch(permz,eachsubjz,pulltopchan,2) = fulldatasetloso(eachsubjz,IE(pulltopchan),2);
                sCOI(permz,eachsubjz,pulltopchan) = IE(pulltopchan);
            end
        end
        %variable of extracted data for use in statistical anlayses for LOSO method
        leftoutsubpull(permz,eachsubjz,1)= nanmean(leftoutsubpullch(permz,eachsubjz,:,1),3); %average data across all channels of interest
        leftoutsubpull(permz,eachsubjz,2)= nanmean(leftoutsubpullch(permz,eachsubjz,:,2),3);
        %clear held out sub vars
        clear meanleftoutsubdiff meanleftoutsub leaveonsubout fulldatasetloso ME IE
    end
end

%%
%%run-level data for within-subjects fcoi methods
iCOIx = zeros(permis,numbsubs,topchan,8);  %set up variable for COIs for leave-one-run-out
for permzz = 1:permis
    for eachsubzz = 1:numbsubs
        insubdata = subrundata(sublist==subnumbs(eachsubzz),:,:);%isolate all runs for a given subject(runs x channels x conds)
        whichsub(permzz,eachsubzz) = subnumbs(eachsubzz);%find subject number
        numbrruns(eachsubzz) = numel(insubdata(:,1,1)); %find number of runs for this subject
        insubpermdata=insubdata;%change variable name for use below
        oddeven = mod(1:numbrruns(eachsubzz),2);%finding odd/even division of runs using mod
        %LEAVE-ONE-RUN-OUT method
        if eachsubzz==1
            startline(eachsubzz) =1;%find index of staring run for this subject in full dataset
        elseif eachsubzz > 1
            startline(eachsubzz) = (sum(numbrruns(1:eachsubzz))- (numbrruns(eachsubzz)))+1;
        end
        for iii = 1:numel(insubpermdata(:,1,1)) % for 1 through # of runs for this subject
            if numbrruns(eachsubzz) > 1 % if more than 1
                loossubdata = insubpermdata;%make temp. dataset for holding out.
                loossubdata(iii,:,:) = [];%delete run of interest from temp. dataset
                mloorunz = nanmean(loossubdata,1);%find mean for each channel for all other runs for each%condidion
                %find differences between conditions
                if activeselect==1   %if focus is active condition
                    mloorunzdiff = mloorunz(1,:,1)-mloorunz(1,:,2);
                elseif activeselect==2%if focus is control condition
                    mloorunzdiff = mloorunz(1,:,2)-mloorunz(1,:,1);
                end
                %find channel that is most responsive to target condition(i.e., active)by sorting diff values
                [MSS, IZZ] = sort(mloorunzdiff,'descend','MissingPlacement','last'); %MSS=diff values; IZZ=channels
                %%extract data from run left out using channel identified in rest of dataset
                for pulltopchaneachrun = 1:topchan
                    if isnan(MSS(1,pulltopchaneachrun))==1 %if top diff is a NaN--there is a problem
                        leftoutdatapullch(iii,pulltopchaneachrun,1)= NaN;
                        leftoutdatapullch(iii,pulltopchaneachrun,2)= NaN;
                        iCOIx(permzz,eachsubzz,pulltopchaneachrun,iii) = NaN;
                    elseif ~isnan(MSS(1,pulltopchaneachrun))==1% if top diff is actua value/not a NaN
                        leftoutdatapullch(iii,pulltopchaneachrun,1) = insubpermdata(iii,IZZ(pulltopchaneachrun),1);%pull all channel data for active cond.
                        leftoutdatapullch(iii,pulltopchaneachrun,2) = insubpermdata(iii,IZZ(pulltopchaneachrun),2);%pull all channel data for control cond.
                        %if pulled channels for both conditions are NaN, then the COI is NaN; implemented for consistency stats
                        if isnan(leftoutdatapullch(iii,pulltopchaneachrun,1))==1 && isnan(leftoutdatapullch(iii,pulltopchaneachrun,2))==1%
                            iCOIx(permzz,eachsubzz,pulltopchaneachrun,iii) = NaN;
                        else
                            iCOIx(permzz,eachsubzz,pulltopchaneachrun,iii) = IZZ(pulltopchaneachrun);
                        end
                    end
                end
                %averaging if there are multiple channels in fCOI def/does nothing if not
                leftoutdatapull(iii,1) = nanmean(leftoutdatapullch(iii,:,1),2);
                leftoutdatapull(iii,2) = nanmean(leftoutdatapullch(iii,:,2),2);
            elseif numbrruns(eachsubzz) < 2 %if subject has less than 2 runs, this subject can't be used
                leftoutdatapull(iii,:) = NaN;%NaN out extraction
            end
            %clear temp. diff value variables for next loop
            clear mloorunz mloorunzdiff IZZ MSS
        end
        %variable of extracted values for LORO method for statistical analysis
        indsubdataCOI(permzz,eachsubzz,1) = nanmean(leftoutdatapull(:,1));%active condition%average of all indepdently-defined data from all runs for this subject
        indsubdataCOI(permzz,eachsubzz,2) = nanmean(leftoutdatapull(:,2));%control condition
        %clear temp. datasets with held out data  for next loop
        clear leftoutdatapull leftoutdatapullch loossubdata
        %
        %SPLIT-HALF method
        for halves = 1:2 %for each half
            if halves ==1 %half1
                halfsubdata = squeeze(nanmean(insubpermdata(oddeven==0,:,:),1));%half1 use even runs to define COI
                Hdiff = halfsubdata(:,1)-halfsubdata(:,2);%calculate differences between conditions
                [HM, HI] = sort(Hdiff,'descend','MissingPlacement','last');%sort differences and corresponding channels
                if numel(topchanrange)==1 %for single channel methods
                    halfpull(permzz,eachsubzz,halves,:) = squeeze(nanmean(insubpermdata(oddeven==1,HI(topchanrange),:)));%pull data from odd runs using COI identified above
                    shCOI(permzz,eachsubzz,topchanrange,1) = HI(topchanrange);%log COI used
                elseif numel(topchanrange)>1 %for multichannel methods
                    halfpull(permzz,eachsubzz,halves,:) = mean(squeeze(nanmean(insubpermdata(oddeven==1,HI(topchanrange),:))));
                    shCOI(permzz,eachsubzz,topchanrange,1) = HI(topchanrange);%log COIs used
                end
            elseif halves ==2 %half2
                halfsubdata = squeeze(nanmean(insubpermdata(oddeven==1,:,:),1));%half2 use odd runs to define COI
                Hdiff = halfsubdata(:,1)-halfsubdata(:,2);%calculate differences between conditions
                [HM, HI] = sort(Hdiff,'descend','MissingPlacement','last');%sort differences and corresponding channels
                if numel(topchanrange)==1 %for single channel methods
                    halfpull(permzz,eachsubzz,halves,:) = squeeze(nanmean(insubpermdata(oddeven==0,HI(topchanrange),:))); %pull data from even runs using COI identified above
                    shCOI(permzz,eachsubzz,topchanrange,2) = HI(topchanrange);%log COI used
                elseif numel(topchanrange)>1 %for multichannel methods
                    halfpull(permzz,eachsubzz,halves,:) = mean(squeeze(nanmean(insubpermdata(oddeven==0,HI(topchanrange),:))));
                    shCOI(permzz,eachsubzz,topchanrange,2) = HI(topchanrange);%log COIs used
                end
            end
            %clear temp variables
            clear halfsubdata Hdiff HM  HI
        end
        %create variable of extracted values using split-half mehtod for further analaysis
        datahalves(permzz,eachsubzz,:) = squeeze(mean(halfpull(permzz,eachsubzz,:,:),3)); %averaging each halves per sub.
        %
        %ITERATIVE CONTRASTS METHOD
        half1sub=squeeze(nanmean(insubpermdata(oddeven==0,:,:),1));%define half1 as even runs
        half2sub=squeeze(nanmean(insubpermdata(oddeven==1,:,:),1));%define half2 as odd runs
        %define 4 unique contrasts
        C(1,:) = half1sub(:,1)-half1sub(:,2); %active half1 vs. control half1
        C(2,:) = half1sub(:,1)-half2sub(:,2); %active half1 vs. control half2
        C(3,:) = half2sub(:,1)-half2sub(:,2); %active half2 vs. control half2
        C(4,:) = half2sub(:,1)-half1sub(:,2); %active half2 vs. control half1
        for contrastz= 1:4 %for each unique contrast, find differences and sort to find/use as COI definition for left out data
            [HMM, HII] = sort(C(contrastz,:),'descend','MissingPlacement','last');%sort diff in contrast
            ICCOI(permzz,eachsubzz,contrastz,:)=HII(topchanrange); %record top COI(s) for each constrast.
            %create vector for logical channel index
            indy=1:10==0;
            for eachel = 1:topchan
                indy(HII(eachel))=1;%make channel index of top chan(s) to pull from
            end
            %extract data from remaining/held out data for each contrast using COI(s)identified above
            if contrastz==1
                icdata(permzz,eachsubzz,contrastz,1) = nanmean(half2sub(indy==1,1));%pull active data for half2
                icdata(permzz,eachsubzz,contrastz,2) = nanmean(half2sub(indy==1,2));%pull control data for half2
            elseif contrastz==2
                icdata(permzz,eachsubzz,contrastz,1) = nanmean(half2sub(indy==1,1));%pull active data for half2
                icdata(permzz,eachsubzz,contrastz,2) = nanmean(half1sub(indy==1,2));%pull control data for half1
            elseif contrastz==3
                icdata(permzz,eachsubzz,contrastz,1) = nanmean(half1sub(indy==1,1));%pull active data for  half1
                icdata(permzz,eachsubzz,contrastz,2) = nanmean(half1sub(indy==1,2));%pull control data for half1
            elseif contrastz==4
                icdata(permzz,eachsubzz,contrastz,1) = nanmean(half1sub(indy==1,1));%pull active data for  half1
                icdata(permzz,eachsubzz,contrastz,2) = nanmean(half2sub(indy==1,2));%pull control data for half2
            end
            indy=1:10==0;%reset channel index
        end
        %create variable of extracted values for iterative contrasts analysis
        %average pulled data for each subject/across all COIs
        Micdata(permzz,eachsubzz,:) = squeeze(nanmean(icdata(permzz,eachsubzz,:,:),3));%average data from extracted contrasts for this subject
        %clear relevant vars
        clear half1sub half2sub C HMM HII oddeven
    end
end
iCOIx(iCOIx==0)=NaN;  %changes remaining zeros in iCOIx matrix to NaN for later consistency metrics

%%
%inferential statistical analyses:
%FIXED Array method
%analyze active vs. control at each channel
channels = 1:numel(subchandata(1,:,1));
channels = channels';
for eachchannel = 1:numel(channels)%for each channel
    %run ttest
    [h,p,ci,stats]=ttest(subchandata(:,eachchannel,1),subchandata(:,eachchannel,2),'Tail','right');%one-tailed for active>control
    chanp(eachchannel,:) = p;%extract p
    chant(eachchannel,:) = stats.tstat;%extract t
    chanCohensD(eachchannel) = nanmean(subchandata(:,eachchannel,1)- subchandata(:,eachchannel,2)) /...
        nanstd(subchandata(:,eachchannel,1)- subchandata(:,eachchannel,2));%calcuate Cohen's D/effect size
end
findsigchans = chanp < .05;%find significant channels without correction, 10 channels 
findsigchanscorrected = chanp < .005;%find signficiant channels after bonferroni correction, 10 channels

%fcoi methods stats
statsperm = 1; %variable formatting for use with permuation tests-not needed for these analyses
%stats for LEAVE ONE SUBJECT OUT method
[h,p,ci,stats]=ttest(leftoutsubpull(statsperm,:,1),leftoutsubpull(statsperm,:,2),'Tail','right');%one-tailed, biased to face.
losop(statsperm) = p;
losot(statsperm) = stats.tstat;
losoCohensD(statsperm) = nanmean(leftoutsubpull(statsperm,:,1)-leftoutsubpull(statsperm,:,2)) /...
    nanstd(leftoutsubpull(statsperm,:,1)-leftoutsubpull(statsperm,:,2));
%stats for LEAVE ONE RUN OUT method
[h,p,ci,stats] = ttest(indsubdataCOI(statsperm,:,1),indsubdataCOI(statsperm,:,2),'Tail','right');
lorop(statsperm) = p;
lorot(statsperm) = stats.tstat;
loroCohensD(statsperm) = nanmean(indsubdataCOI(statsperm,:,1)- indsubdataCOI(statsperm,:,2)) /nanstd(indsubdataCOI(statsperm,:,1)- indsubdataCOI(statsperm,:,2));
%stats for SPLIT-HALF method
[h,p,ci,stats] = ttest(datahalves(statsperm,:,1),datahalves(statsperm,:,2),'Tail','right');
shp(statsperm) = p;
sht(statsperm) = stats.tstat;
shCohensD(statsperm) = nanmean(datahalves(statsperm,:,1)- datahalves(statsperm,:,2)) /nanstd(datahalves(statsperm,:,1)- datahalves(statsperm,:,2));
%stats for ITERATIVE CONTRASTS method
[h,p,ci,stats] = ttest(Micdata(statsperm,:,1),Micdata(statsperm,:,2),'Tail','right');
pp(statsperm) = p;
pt(statsperm) = stats.tstat;
pCohensD(statsperm) = nanmean(Micdata(statsperm,:,1)- Micdata(statsperm,:,2)) /nanstd(Micdata(statsperm,:,1)- Micdata(statsperm,:,2));


%%
%Consistency metrics/analysis
%within-ss LORO consistency
for eachsubzzz = 1:numel(unique(sublist))%for each subject
    A = squeeze(iCOIx(statsperm,eachsubzzz,:,:));%reduce COI variable down to be workable
    countsofcoi = [min(min(A)):max(max(A));histcounts(A)]';%count the frequency of each channel being defined as COI
    [B,I] = sort(countsofcoi(:,2),'descend');%sort frequencies
    percentagesameICOI(statsperm,eachsubzzz)= sum(countsofcoi(I(topchanrange),2))/sum(~isnan(A),'all');%%proportion of top X most frequ. channels/total COIs identified for subject
    clear A B I countsofcoi
end
MoverallpercentagetopCOI(statsperm,:) = nanmean(percentagesameICOI(statsperm,:),2); %averages proportions across all subjects
%
%within-ss splif-half consistency==number that match / total number
%(assigns 0% to non-matches rather than 50%)
if topchan == 1 %for single channel defs
    MoverallproportionsameSH(statsperm) = sum(shCOI(statsperm,:,1,1)==shCOI(statsperm,:,1,2))/numel(shCOI(statsperm,:,1,1)); %total matches between pairs/total pairs
elseif topchan == 2 %for 2 channel defs
    MoverallproportionsameSH(statsperm) = sum(shCOI(statsperm,:,:,2)==shCOI(statsperm,:,1,1)|shCOI(statsperm,:,:,2)==shCOI(statsperm,:,2,1),'all')/ numel(shCOI(statsperm,:,:,1));%total matches from combination of 2 pairs/total pairs
elseif topchan==3 %for 3 channel defs.; %total matches from combination of 3 pairs/total pairs
    MoverallproportionsameSH(statsperm) = sum(shCOI(statsperm,:,:,2)==shCOI(statsperm,:,1,1)|shCOI(statsperm,:,:,2)==shCOI(statsperm,:,2,1)|shCOI(statsperm,:,:,2)==shCOI(statsperm,:,3,1),'all')/ numel(shCOI(statsperm,:,:,1));
end
%
%within-ss, iterative contrasts consistency
for eachsubzzzz = 1:numel(ICCOI(statsperm,:,1,1))
    AA = squeeze(ICCOI(statsperm,eachsubzzzz,:,:)); %squeezes into contrasts(4) x nunber of cois
    countsofcoII = [min(min(AA)):max(max(AA));histcounts(AA)]';%frequency each channel chosen as COI
    [B,I] = sort(countsofcoII(:,2),'descend');%sort
    topchanz(:,eachsubzzzz) = countsofcoII(I(topchanrange));%identify top X channels
    percentagesameIC(statsperm,eachsubzzzz)=sum(countsofcoII(I(topchanrange),2))/numel(AA);%proportion that match the top X channels/total COIs chosen
    clear AA B I
end
MoverallpercentagetopIC(statsperm,:) = nanmean(percentagesameIC(statsperm,:));%
%
%between-subjects LOSO consistency
countzofcoi= [min(min(sCOI(statsperm,:,:))):max(max(sCOI(statsperm,:,:)));histcounts(sCOI(statsperm,:,:))]';%frequency of times each channel selected as COI
[B,I] = sort(countzofcoi(:,2),'descend');%sort
topchanz = countzofcoi( I(topchanrange));%identify the top channels
if topchan==1 %if single channel def.
    MoverallproportionsameLOSO(statsperm) =sum(sCOI(statsperm,:)==topchanz(1), 'all')/numel(sCOI(statsperm,:));%proportion of matches to top channels/total COIs chosen
elseif topchan==2%for 2 channel defs.
    MoverallproportionsameLOSO(statsperm) = sum(sCOI(statsperm,:)==topchanz(1)|sCOI(statsperm,:)==topchanz(2),'all')/numel(sCOI(statsperm,:));%sum that match top X/total
elseif topchan==3%for 3 channel defs.
    MoverallproportionsameLOSO(statsperm) = sum(sCOI(statsperm,:)==topchanz(1)|sCOI(statsperm,:)==topchanz(2)|sCOI(statsperm,:)==topchanz(3),'all')/numel(sCOI(statsperm,:));%sum that match top X/total
end
clear countzofcoi topchanz

%displaying main statistical results as output in command window%
display('fixed array results')
chant
chanp
chanCohensD'
display('LOSO results')
losot
losop
losoCohensD
display('LORO results')
lorot
lorop
loroCohensD
display('split-half results')
sht
shp
shCohensD
display('iterative contrasts results')
pt
pp
pCohensD
