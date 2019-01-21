%% INITIAL CONDITIONS & PARAMETRIZATION %%
clear all;close all;clc;
for rateOfParentalConfusion = 0:0.1:1 % Rate of parental confusion among males.
    clearvars -except rateOfParentalConfusion;close all;clc;
    numOfMales           = 100; % Initial number of males in the population
    numOfFemales         = 100; % Initial number of females in the population
    meanAgeMales         = 15*360; % mean for the age distribution of males (days)
    meanAgeFemales       = 15*360; % mean for the age distribution of females (days)
    lifeExpectancy       = 33*360; % Life expectancy if 33 years
    stdAgeMales          = 10*360; % Standart deviation for the age of males
    stdAgeFemales        = 10*360; % Standart deviation for the age of females
    fertilityAgeMales    = 12*360; % Age when males become fertile (days)
    fertilityAgeFemales  = 12*360; % Age when females become fertile (days)
    mensturalCycleLength = 28; % Menstural cycle length of females (days)
    dayOvulationBegins   = 11; % Day that the female ovulation begins (therefore can get pregnant afterwards)
    dayOvulationEnds     = 21; % Day that the female ovulation ends (therefore can't get pregnant afterwards)
    gestationPeriod      = 280; % Period of Gestation in humans (days)
    lactationalAmenorrheaLength    = 6*30; % Period of Postpartum amenorrhea given brestfeeding (lactation) (days)
    nonlactationalAmenorrheaLength = 2*30; % Period of Postpartum amenorrhea given NO brestfeeding (lactation) (days)
    p                              = 1-(rateOfParentalConfusion)^(1/nonlactationalAmenorrheaLength); % Probability of being murdered per day
    totalDaysOfSimualtion          = 360*50; % Total days of simulation
    probabilityOfPregnancy         = 1;% prob that the female will get pregnant after copulation
    % RANDOMIZE THE AGE OF MALES / FEMALES (POISSON DIST) %%
    agesOfMales   = abs(normrnd(meanAgeMales,stdAgeMales,[1 numOfMales]));
    agesOfFemales = abs(normrnd(meanAgeFemales,stdAgeFemales,[1 numOfFemales]));
    %% RANDOMIZE THE DAY OF MENSTURAL CYCLE THAT FEMALES ARE AT (UNIFORM DIST) %%
    dayAtMensturalCycle = randi(28,[1 numOfFemales]); % Females ovulate only between the days 11-21, and if the female is below 12 years old, no matter at which day she is, she won't copulate / get pregnant
    statusPregnancy     = zeros(1,numOfFemales);
    statusLAM           = zeros(1,numOfFemales);
    dayPregnancy        = zeros(1,numOfFemales);
    dayLAM              = zeros(1,numOfFemales);
    mumIndxMales        = zeros(1,numOfMales); % initially we don't know the mothers
    mumIndxFemales      = zeros(1,numOfFemales); % initially we don't know the mothers
    maleMatrix   = [(1:numOfMales)' agesOfMales' mumIndxMales'];
    femaleMatrix = [(1:numOfFemales)' agesOfFemales' dayAtMensturalCycle' statusPregnancy' statusLAM' dayPregnancy' dayLAM' mumIndxFemales'];
    %%
    lastMaleIdx   = numOfMales;
    lastFemaleIdx = numOfFemales;
    killedkids    = [];
    for day=1:totalDaysOfSimualtion
        day
        %% CHECK THE ALREADY PREGNANT FEMALES
        alreadyPregnantFemaleIndexes =  femaleMatrix(femaleMatrix(:,4)==1,1);
        for f=1:length(alreadyPregnantFemaleIndexes)
            alreadyPregnantFemaleIdx = alreadyPregnantFemaleIndexes(f);
            if(femaleMatrix(femaleMatrix(:,1)==alreadyPregnantFemaleIdx,6)==gestationPeriod) % GIVING BIRTH -> PUT THE FEMALE IN LAM, ADD THE KID TO F OR M
                genderOfNewBorn = randsample([0,1],1); % (female = 0, male = 1) ->50% chance, can be changed
                if(genderOfNewBorn==0) % it's a female
                    femaleMatrix  = [femaleMatrix; lastFemaleIdx+1 0 1 0 0 0 0 alreadyPregnantFemaleIdx];
                    lastFemaleIdx = lastFemaleIdx+1;
                else % it's a male
                    maleMatrix   = [maleMatrix; lastMaleIdx+1 0 alreadyPregnantFemaleIdx];
                    lastMaleIdx  = lastMaleIdx+1;
                end
                femaleMatrix(femaleMatrix(:,1)==alreadyPregnantFemaleIdx,3)=0;% SET THE DAY AT MENSTURAL CYCLE TO 0 (JUST ENTERED LAM)
                femaleMatrix(femaleMatrix(:,1)==alreadyPregnantFemaleIdx,4)=0;% SET PREG. STATUS TO 0
                femaleMatrix(femaleMatrix(:,1)==alreadyPregnantFemaleIdx,5)=1;% SET LAM STATUS TO 1
                femaleMatrix(femaleMatrix(:,1)==alreadyPregnantFemaleIdx,6)=0;% SET PREG. DAYS TO 0
                femaleMatrix(femaleMatrix(:,1)==alreadyPregnantFemaleIdx,7)=1;% SET LAM DAYS TO 1
            else % IF NOT GIVING BIRTH, INCREASE THE PREGNANCY DAY BY ONE
                femaleMatrix(femaleMatrix(:,1)==alreadyPregnantFemaleIdx,6)=femaleMatrix(femaleMatrix(:,1)==alreadyPregnantFemaleIdx,6)+1;
            end
        end
        %% CHECK THE ALREADY LAM FEMALES
        alreadyLAMFemaleIndexes =  femaleMatrix(femaleMatrix(:,5)==1,1);
        for f=1:length(alreadyLAMFemaleIndexes)
            alreadyLAMFemaleIdx = alreadyLAMFemaleIndexes(f);
            if(femaleMatrix(femaleMatrix(:,1)==alreadyLAMFemaleIdx,7)==lactationalAmenorrheaLength) % GETTING OUT OF LAM -> SET THE DAY OF MENS. CYCLYE TO 1, SET THEIR STATUS AS NEITHER PREGNANT NOR LAM
                femaleMatrix(femaleMatrix(:,1)==alreadyLAMFemaleIdx,3)=1;% SET THE DAY AT MENSTURAL CYCLE TO 1 (JUST GOT OUT OF LAM)
                femaleMatrix(femaleMatrix(:,1)==alreadyLAMFemaleIdx,7)=0;% SET THE DAY AT LAM TO 0 (JUST GOT OUT OF LAM)
                femaleMatrix(femaleMatrix(:,1)==alreadyLAMFemaleIdx,5)=0;% SET THE LAM STATUS TO 0
            else % IF NOT GETTING OUT OF LAM
                attacked = randsample([0,1],1,true,[1-p p]);
                if(attacked==1) % IF ATTACKED
                    % THE KID OF THIS MOM GETS KILLED (randomly selected, but has to be younger then 6 mo)
                    femaleKidsOfMum = femaleMatrix((femaleMatrix(:,8)==alreadyLAMFemaleIdx)&(femaleMatrix(:,2)<lactationalAmenorrheaLength),1);
                    maleKidsOfMum   = maleMatrix((maleMatrix(:,3)==alreadyLAMFemaleIdx)&(maleMatrix(:,2)<lactationalAmenorrheaLength),1);
                    if(~isempty(femaleKidsOfMum) && ~isempty(maleKidsOfMum)) %if the mum has both male and female kids
                        gender2beKilled = randsample([0,1],1); % (female = 0, male = 1) ->50% chance, can be changed
                        if(gender2beKilled==0)
                            female2bekilled = randsample([femaleKidsOfMum femaleKidsOfMum],1); %technical, need to double to make it seem like a vector in case it's only one kid
                            femaleMatrix(femaleMatrix(:,1)==female2bekilled,:)=[]; % KILLED :(
                            killedkids = [killedkids; day 0 female2bekilled];
                        else
                            male2bekilled = randsample([maleKidsOfMum maleKidsOfMum],1);
                            maleMatrix(maleMatrix(:,1)==male2bekilled,:)=[]; % KILLED :(
                            killedkids = [killedkids; day 1 male2bekilled];
                        end
                    elseif(~isempty(femaleKidsOfMum) && isempty(maleKidsOfMum))
                        female2bekilled = randsample([femaleKidsOfMum femaleKidsOfMum],1);
                        femaleMatrix(femaleMatrix(:,1)==female2bekilled,:)=[]; % KILLED :(
                        killedkids = [killedkids; day 0 female2bekilled];
                    elseif(isempty(femaleKidsOfMum) && ~isempty(maleKidsOfMum))
                        male2bekilled = randsample([maleKidsOfMum maleKidsOfMum],1);
                        maleMatrix(maleMatrix(:,1)==male2bekilled,:)=[]; % KILLED :(
                        killedkids = [killedkids; day 1 male2bekilled];
                    end
                    % FEMALE GETS OUT OF LAM PREMATURELY, SET THE LAM DAYS
                    % (lactationalAmenorrheaLength-nonlactationalAmenorrheaLength)
                    femaleMatrix(femaleMatrix(:,1)==alreadyLAMFemaleIdx,3)=1;% SET THE DAY AT MENSTURAL CYCLE TO 1 (JUST GOT OUT OF LAM)
                    femaleMatrix(femaleMatrix(:,1)==alreadyLAMFemaleIdx,7)=lactationalAmenorrheaLength-nonlactationalAmenorrheaLength;% INCREASE THE DAY AT LAM TO 0 (JUST GOT OUT OF LAM)
                    femaleMatrix(femaleMatrix(:,1)==alreadyLAMFemaleIdx,5)=0;% SET THE LAM STATUS TO 0
                else
                    % IF NOT ATTACKED INCREASE LAM DAYS BY 1
                    femaleMatrix(femaleMatrix(:,1)==alreadyLAMFemaleIdx,7)=femaleMatrix(femaleMatrix(:,1)==alreadyLAMFemaleIdx,7)+1;
                end
            end
        end
        
        %% CHECK THE AVALIABLE MALES AND FEMALES - GET NEW FEMALES PREGNANT
        avaliableMaleIndexes   = maleMatrix(maleMatrix(:,2)>fertilityAgeMales,1); % indexes of the males avaliable for copulation
        avaliableFemaleIndexes = femaleMatrix((femaleMatrix(:,2)>fertilityAgeFemales)&(femaleMatrix(:,3)>dayOvulationBegins)&(femaleMatrix(:,3)<dayOvulationEnds)&(femaleMatrix(:,4)==0)&(femaleMatrix(:,5)==0),1); % indexes of the females avaliable for pregnancy
        newPregnantFemaleIndexes  = randsample(avaliableFemaleIndexes,round(min(length(avaliableMaleIndexes),length(avaliableFemaleIndexes))*probabilityOfPregnancy));
        for f=1:length(newPregnantFemaleIndexes)
            newlyPregnantFemaleIdx = newPregnantFemaleIndexes(f);
            femaleMatrix(femaleMatrix(:,1)==newlyPregnantFemaleIdx,4)=1; %set the status of being pregnant to 1;
            femaleMatrix(femaleMatrix(:,1)==newlyPregnantFemaleIdx,6)=1; %set the day of pregnancy to 1;
            femaleMatrix(femaleMatrix(:,1)==newlyPregnantFemaleIdx,3)=0; %set the day of menstural cycle to 0;
        end
        %% INCREASE EVERYBODYS AGE BY 1 DAY
        femaleMatrix(:,2) = femaleMatrix(:,2)+1;
        maleMatrix(:,2)   = maleMatrix(:,2)+1;
        %% KILL THE ONES THAT ARE OVER 33 YEARS OLD (random death can also be included)
        femaleMatrix(femaleMatrix(:,2)>lifeExpectancy,:)=[];
        maleMatrix(maleMatrix(:,2)>lifeExpectancy,:)=[];
        %% INCREASE THE DAY AT MENSTURAL CYCLE FOR NON-PREGNANT NON-LAM FEMALES
        femaleMatrix((femaleMatrix(:,4)==0)&(femaleMatrix(:,5)==0),3)=femaleMatrix((femaleMatrix(:,4)==0)&(femaleMatrix(:,5)==0),3)+1;
        femaleMatrix(:,3) = mod(femaleMatrix(:,3),29); %set the menstural cycle day to modula 29
        keepSize(day,:)   = [size(femaleMatrix,1) size(maleMatrix,1)];
    end
    save(['COEVOLUTION_DATA_' num2str(10*rateOfParentalConfusion)]) %save the dataset 
end













