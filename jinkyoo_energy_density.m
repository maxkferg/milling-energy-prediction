%clear all


process_ID=[1];

%base covariance function scale parameters
Hyp = [log(1000),log(6000),log(6),log(3),log(3),log(3),log(3),log(3),log(3),log(3),log(1)];


%cut with feed
Feature_set{1} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{2} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{3} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{4} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{5} = [1 2 3   4 5 6 7   8 9 10 ];
Feature_set{6} = [1 2 3   4 5 6 7   8 9 10 ];

%air cut
Feature_set{7} = [1 2 3   4 5 6 7   8 9 10 ];

%Rapid motion
Feature_set{8} = [1 2 3   4 5 6 7   8 9 10 ];

%Dwell
Feature_set{9} = [1 2 3   4 5 6 7   8 9 10 ];

%Etc
Feature_set{10} = [1 2 3   4 5 6 7   8 9 10 ];




load('data/Training1.mat')
load('data/Training2.mat')
load('data/Training3.mat')
load('data/Training4.mat')
load('data/Training5.mat')
load('data/Training6.mat')
load('data/Training7.mat')
load('data/Training8.mat')
load('data/Training9.mat')
load('data/Training10.mat')
load('data/Training11.mat')
load('data/Training12.mat')
load('data/Training13.mat')
load('data/Training14.mat')
load('data/Training15.mat')
load('data/Training16.mat')
load('data/Training17.mat')
load('data/Training18.mat')

D = [Training1;Training2;Training3;Training4;Training5;Training6;Training7;Training8;Training9;Training10;Training11;Training12;Training13;Training14;Training15;Training16;Training17;Training18];
%D = [Training1;Training2;Training3;Training4;Training5;Training6;Training7;Training8;Training9];
%D = [Training10;Training11;Training12;Training13;Training14;Training15;Training16;Training17;Training18];



%% extract the field
energy = cell2mat(D(:,9)); %energy consumption
duration = cell2mat(D(:,10)); %duration of operation
feed = cell2mat(D(:,11)); %duration of operation
spindle = cell2mat(D(:,12)); %spindle speed
length_cut_X = abs(cell2mat(D(:,19))); %code dx
length_cut_Y = abs(cell2mat(D(:,20))); %code dy
length_cut_Z = abs(cell2mat(D(:,23))); %code dy
length_cut_XY = cell2mat(D(:,21)); %code length_cut
length_cut_XYZ = cell2mat(D(:,31)); %code length_cut
actual_dx = cell2mat(D(:,25)); %actual dx
actual_dy = cell2mat(D(:,26)); %actual dy
actual_length_cut = sqrt(actual_dx.^2+actual_dy.^2);
depth_cut = cell2mat(D(:,27)); %Depth of cut
area_cut = cell2mat(D(:,28)); %Depth of cut
volume_cut = cell2mat(D(:,29)); %Depth of cut

%ratio cut
ratio_cut = actual_length_cut./length_cut_XY;
for i=1:length(ratio_cut)
    if ratio_cut(i)>0;
        ratio_cut(i)=ratio_cut(i);
    else
        ratio_cut(i)=0;
    end
end

%cut_method=zeros(length(y),1);
for i=1:length(energy)
    
    if strcmp(D{i,24},'Conventional')
        cut_method(i,:) = [1,0,0];
    elseif strcmp(D{i,24},'Climb')
        cut_method(i,:) = [0,1,0];
    else
        cut_method(i,:)=[0,0,1];
    end
end

%cut_direction=zeros(length(y),1);
for i=1:length(energy)
    if abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))==0
        cut_direction(i,:) = [1,0,0,0];
    elseif abs(cell2mat(D(i,19)))==0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,:) = [0,1,0,0];
    elseif abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,:)=[0,0,1,0];
    else
        cut_direction(i,:)=[0,0,0,1];
    end
    
end


%opeartion
for i=1:length(energy)
    
    if strcmp(D{i,32},'Face Milling')
        type_operation(i,1) = 1;
    elseif strcmp(D{i,32},'Contouring')
        type_operation(i,1) = 2;
    elseif strcmp(D{i,32},'Slotting')
        type_operation(i,1) = 3;
    elseif strcmp(D{i,32},'Pocketing')
        type_operation(i,1) = 4;
    else
        type_operation(i,1) = 0;
    end
    
end


%label
for i=1:length(energy)
    if strcmp(D{i,30},'Cut with Feed')
        label(i,1) = 1;
    elseif strcmp(D{i,30},'Plunge with feed')
        label(i,1) = 2;
        
        
        
    elseif strcmp(D{i,30},'Air-Cut')
        label(i,1) = 3;
    elseif strcmp(D{i,30},'Air-Cut in Z while plunging')
        label(i,1) = 4;
    elseif strcmp(D{i,30},'Air-cut in Z while retracting')
        label(i,1) = 5;
        
        
        
    elseif strcmp(D{i,30},'No Cut - Rapid motion')
        label(i,1) = 6;
        
    elseif strcmp(D{i,30},'Dwell')
        label(i,1) = 7;
        
    else
        label(i,1) = 0;
    end
end






for i=1:length(energy)
    % cut with feed
    if (label(i) == 1) %(cut in x-y direction)
        
        if (type_operation(i) == 1)     % face milling
            ID(i,1) = 1;
        elseif (type_operation(i) == 2) % contouring
            ID(i,1) = 2;
        elseif (type_operation(i) == 3) % splitting
            ID(i,1) = 3;
        elseif (type_operation(i) == 4) % pocketing
            ID(i,1) = 4;
        else
            ID(i,1) = 0; %
        end
    elseif (label(i)==2)
        ID(i,1) = 5;   %Plunge with feed
    elseif (label(i) == 3  || label(i) == 4 || label(i) == 5) % air cut
        ID(i,1) = 6;
    elseif (label(i) == 6) %rapid motion
        ID(i,1) = 7;
    elseif (label(i) == 7) %dwell
        ID(i,1) = 8;
    else
        ID(i,1) = 9; %no-labeling
    end
end






input= [feed,...,        %1
    spindle,...,        %2
    depth_cut,...,      %3
    cut_direction,...,  %4 5 6 7
    cut_method,...,     %8 9 10
    length_cut_XYZ,..., %11
    ID,...,             %12
    duration];          %13
output =energy;
%density=energy./length_cut_XYZ;

%total data set
X = input;
E = energy;
for i=1:length(E)
    if ID(i) == 8 %dwell
        L(i,1)=1;
    else
        L(i,1)=X(i,11);
    end
end
Y = E./L;



index_zero = find(Y==0);
index_zero=union([index_zero], [index_zero+1])
index_zero=union([index_zero], [index_zero-1])

%clean up data
clean_up_index = setdiff([1:length(energy)],[index_zero]);
X = X(clean_up_index,:);
ID = ID(clean_up_index);
L = L(clean_up_index);
E = E(clean_up_index);
Y = Y(clean_up_index);

%clean up data
clean_up_index = find(Y > 0 & Y < 3000 & L > 0);
X = X(clean_up_index,:);
ID = ID(clean_up_index);
L = L(clean_up_index);
E = E(clean_up_index);
Y = Y(clean_up_index);



%training individul prediction function
for I=1:length(process_ID)
    
    feature_index = Feature_set{I}
    process_id = process_ID(I);
    
    %Select data depending on the process ID
    if process_id == 7 %rapid motion
        index_job{I} = find(X(:,12)==process_id );
    else                                             %time           %non-zero cut
        index_job{I} = find(X(:,12)==process_id & X(:,13)>2 & X(:,11)>0); %cutting related
    end
    
    %selection data
    X_training=X(index_job{I},feature_index);
    E_training=E(index_job{I});
    L_training=L(index_job{I});
    T_training=X(index_job{I},13);
    Y_training=E_training./L_training;
    
    
    
    %Learning GP
    hyp.cov = [Hyp(feature_index),0];
    meanfunc = {@meanSum, {@meanLinear, @meanConst}}; %hyp.mean = zeros(size(X,2)+1,1);
    likfunc = @likGauss; sn = 500; hyp.lik = log(sn);
    covfunc = @covSEard;
    % Comment on subsequent runs
    %hyp2 = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, X_training, Y_training)
    f{I}=hyp2
    STD{I}=std(X_training)
end








inputFiles_blind{1}='data/blind_test_accurate.mat'
inputFiles_blind{2}='data/blind_test_intermediate.mat'
inputFiles_blind{3}='data/blind_test_bad.mat'

for I=1:length(inputFiles_blind)
    
    %prediction part
    clear cut_method cut_direction type_operation label ID input E_predict S_predict Y_predict L_blind Y_blind X_blind
    load(inputFiles_blind{I})
    
    D = [Data];
    
    
    
    %extract fetures
    energy = cell2mat(D(:,9)); %answer
    duration = cell2mat(D(:,10)); %duration
    feed = cell2mat(D(:,11)); %duration of operation
    spindle = cell2mat(D(:,12)); %spindle speed
    length_cut_X = abs(cell2mat(D(:,19))); %code dx
    length_cut_Y = abs(cell2mat(D(:,20))); %code dy
    length_cut_Z = abs(cell2mat(D(:,23))); %code dy
    length_cut_XY = cell2mat(D(:,21)); %code length_cut
    length_cut_XYZ = cell2mat(D(:,31)); %code length_cut
    actual_dx = cell2mat(D(:,25)); %actual dx
    actual_dy = cell2mat(D(:,26)); %actual dy
    actual_length_cut = sqrt(actual_dx.^2+actual_dy.^2);
    depth_cut = cell2mat(D(:,27)); %Depth of cut
    volume_cut = cell2mat(D(:,29)); %Depth of cut
    
    ratio_cut = actual_length_cut./length_cut_XY;
    for i=1:length(ratio_cut)
        if ratio_cut(i)>0;
            ratio_cut(i)=ratio_cut(i);
        else
            ratio_cut(i)=0;
        end
    end
    
    %cut_method=zeros(length(y),1);
    for i=1:length(energy)
        
        if strcmp(D{i,24},'Conventional')
            cut_method(i,:) = [1,0,0];
        elseif strcmp(D{i,24},'Climb')
            cut_method(i,:) = [0,1,0];
        else
            cut_method(i,:)=[0,0,1];
        end
    end
    
    %cut_direction=zeros(length(y),1);
    for i=1:length(energy)
        if abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))==0
            cut_direction(i,:) = [1,0,0,0];
        elseif abs(cell2mat(D(i,19)))==0 & abs(cell2mat(D(i,20)))>0
            cut_direction(i,:) = [0,1,0,0];
        elseif abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))>0
            cut_direction(i,:)=[0,0,1,0];
        else
            cut_direction(i,:)=[0,0,0,1];
        end
        
    end
    
    
    %opeartion
    for i=1:length(feed)
        
        if strcmp(D{i,33},'Face Milling')
            type_operation(i,1) = 1;
        elseif strcmp(D{i,33},'Contouring')
            type_operation(i,1) = 2;
        elseif strcmp(D{i,33},'Slotting')
            type_operation(i,1) = 3;
        elseif strcmp(D{i,33},'Pocketing')
            type_operation(i,1) = 4;
        else
            type_operation(i,1) = 0;
        end
    end
    
    %label
    for i=1:length(feed)
        
        if strcmp(D{i,30},'Cut with Feed')
            label(i,1) = 1;
        elseif strcmp(D{i,30},'Plunge with feed')
            label(i,1) = 2;
            
            
            
        elseif strcmp(D{i,30},'Air-Cut')
            label(i,1) = 3;
        elseif strcmp(D{i,30},'Air-Cut in Z while plunging')
            label(i,1) = 4;
        elseif strcmp(D{i,30},'Air-cut in Z while retracting')
            label(i,1) = 5;
            
            
        elseif strcmp(D{i,30},'No Cut - Rapid motion')
            label(i,1) = 6;
        elseif strcmp(D{i,30},'Rapid retract')
            label(i,1) = 6;
            
        elseif strcmp(D{i,30},'Dwell')
            label(i,1) = 7;
        else
            label(i,1) = 0;
        end
    end
    
    
    
    
    
    for i=1:length(energy)
        % cut with feed
        if (label(i) == 1) %(cut in x-y direction)
            
            if (type_operation(i) == 1)     % face milling
                ID(i,1) = 1;
            elseif (type_operation(i) == 2) % contouring
                ID(i,1) = 2;
            elseif (type_operation(i) == 3) % splitting
                ID(i,1) = 3;
            elseif (type_operation(i) == 4) % pocketing
                ID(i,1) = 4;
            else
                ID(i,1) = 0; %
            end
        elseif (label(i)==2)
            ID(i,1) = 5;   %Plunge with feed
        elseif (label(i) == 3  || label(i) == 4 || label(i) == 5) % air cut
            ID(i,1) = 6;
        elseif (label(i) == 6) %rapid motion
            ID(i,1) = 7;
        elseif (label(i) == 7) %dwell
            ID(i,1) = 8;
        else
            ID(i,1) = 9; %no-labeling
        end
    end
    
    
    
    input= [feed,...,        %1
        spindle,...,        %2
        depth_cut,...,      %3
        cut_direction,...,  %4 5 6 7
        cut_method,...,     %8 9 10
        length_cut_XYZ,..., %11
        ID,...,             %12
        duration];          %13
    
    output =energy;
    density=energy./length_cut_XYZ;
    
    
    %total data set
    X_blind = input;
    ID_blind = input(:,12);
    E_blind = energy;
    for i=1:length(E_blind)
        if ID(i) == 8 %dwell
            L_blind(i,1)=1;
        else
            L_blind(i,1)=X_blind(i,11);
        end
    end
    Y_blind = E_blind./L_blind;
    T_blind = input(:,13);
    
    
    
    index_zero = find(Y_blind==0);
    index_zero=union([index_zero], [index_zero+1])
    index_zero=union([index_zero], [index_zero-1])
    
    %clean up data
    clean_up_index = setdiff([1:length(energy)],[index_zero]);
    X_blind = X_blind(clean_up_index,:);
    ID_blind = ID_blind(clean_up_index);
    L_blind = L_blind(clean_up_index);
    E_blind = E_blind(clean_up_index);
    Y_blind = Y_blind(clean_up_index);
    T_blind = T_blind(clean_up_index);
    
    
    
    % %clean up data
    clean_up_index = find(Y_blind>0 & Y_blind < 10000 & L_blind > 0 & ID_blind==1);
    X_blind = X_blind(clean_up_index,:);
    ID_blind = ID_blind(clean_up_index);
    L_blind = L_blind(clean_up_index);
    E_blind = E_blind(clean_up_index);
    Y_blind = Y_blind(clean_up_index);
    T_blind = T_blind(clean_up_index);
    
    
    
    % Preload all of the models from PMML
    % This is a slow operation so we only want to do it once
    % Code block written by Max Ferguson April 14th 2016
    for i=1:length(E_blind)
        process_id = ID_blind(i);
        if process_id<0 || process_id >=10
            continue; % Invalid model number
        elseif exist('models','var') && length(models)>=process_id && ~isempty(models{process_id})
            continue % Already loaded this model
        end
        % Load the correct model file from PMML
        filename = sprintf('models/energy-model-%i.pmml',process_id);
        fprintf('Loading model from %s\n',filename);
        models{process_id} = pmml.GaussianProcess(filename);
    end
    
    
    for i=1:length(E_blind)
        %select the right prediction function
        
        process_id = ID_blind(i);
        
        if (process_id>0 && process_id <9)
            hyp = f{process_id};
            X_training = X(index_job{process_id}, Feature_set{process_id});
            Y_training = Y(index_job{process_id});
            
            [Y_predict(i,1), S_predict(i,1) ] = gp(hyp, @infExact, [], covfunc, likfunc, X_training, Y_training, X_blind(i,Feature_set{process_id}));
            %S_predict(i,1)=S_predict(i,1)+sn;
            E_predict(i,1) = Y_predict(i,1).*L_blind(i);
        else
            Y_predict(i,1)=0;
            S_predict(i,1)=0;
            E_predict(i,1)=0;
            %                 Y_blind(i,1)=0;
            %                 S_blind(i,1)=0;
            %                 E_blind(i,1)=0;
        end
    end
    
    Number_b{I} = length(T_blind);
    Duration_b{I} = mean(T_blind);
    RAE_density_b{I} = mean(abs(Y_predict.*L_blind-Y_blind.*L_blind))/mean(Y_blind.*L_blind);
    %RAE_density_b{I} = mean(abs(Y_predict-Y_blind))/mean(Y_blind);
    predicted_energy_b{I} = sum(Y_predict.*L_blind);
    measured_energy_b{I} = sum(E_blind);
    standard_b{I} = sqrt(sum(S_predict.*L_blind.^2));
    RTE_density_b{I} = (sum(Y_predict.*L_blind)-sum(E_blind))/sum(E_blind);
    
    %display([number,duration,RAE_density*100,predicted_energy/1000,measured_energy/1000,standard,RTE_density*100])
    XX_measured_b{I} = X_blind;
    YY_measured_b{I} = Y_blind;
    EE_measured_b{I} = E_blind;
    LL_measured_b{I} = L_blind;
    
    EE_predict_b{I} = E_predict;
    SS_predict_b{I} = S_predict;
    YY_predict_b{I} = Y_predict;
    
end





















%%no dee to classifiy the spindle speed???????????






%first blind part
XXX = XX_measured_b{1}
YYY = YY_measured_b{1};

feed_range = 0:1:1000;
feed = XXX(:,1);
spindle = XXX(:,2);
depth_cut = XXX(:,3);
cut_direction = XXX(:,4:7);
cut_method = XXX(:,8:10);
density = YYY;

%cut == 1          input(:,4)        input(:,5)         input(:,3)    input(:,2)     input(:,2)    input(:,12)
index_11 = find(cut_direction(:,2)==1 & cut_method(:,1)==1 & depth_cut == 1 & spindle >1400 & spindle <1600 );
index_12 = find(cut_direction(:,2)==1 & cut_method(:,1)==1 & depth_cut == 1 & spindle >2900 & spindle <3100 );
index_13 = find(cut_direction(:,2)==1 & cut_method(:,1)==1 & depth_cut == 1 & spindle >4400 & spindle <4600 );

nrows = length(feed_range);

%feed            spindle                          depth_cut                     cut_direction                   cut_method
XX_11 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),repmat([0 1 0 0],length(feed_range),1), repmat([1 0 0],length(feed_range),1)];
XX_12 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),repmat([0 1 0 0],length(feed_range),1), repmat([1 0 0],length(feed_range),1)];
XX_13 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),repmat([0 1 0 0],length(feed_range),1), repmat([1 0 0],length(feed_range),1)];

% Load the file from a PMML file and plot it on top.
% We choose to load model 1 for facemilling
% We need to convert XX_11, XX_12, XX_13 to the five feature form used by the PMML Model
model = models{1};
XPMML_11 = [XX_11(:,1:3), repmat([2 1],nrows,1)];
XPMML_12 = [XX_12(:,1:3), repmat([2 1],nrows,1)];
XPMML_13 = [XX_13(:,1:3), repmat([2 1],nrows,1)];

%subplot('Position',[0.05 0.1 0.25 0.8]);
figure(1)
subplot('Position',[0.05 0.2 0.8/3 0.75]);
hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_11);
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','b');
plot(feed_range,m,'--b','linewidth',1.5)
% Plot the model trained by Jinkyoo
plot(feed(index_11),density(index_11),'ob','Linewidth',1.5)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11)
plot(feed_range,m,'--r','linewidth',1.5)
legend('Measured','Predicted mean')
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r')
alpha(0.2)
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(150,15,'Spindle speed = 1,500 (RPM)','Interpreter','Latex')
box on


subplot('Position',[0.05*2+0.8/3+0.01 0.2 0.8/3 0.75]);
hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_12);
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','b');
plot(feed_range,m,'--b','linewidth',1.5)
% Plot the model trained by Jinkyoo
plot(feed(index_12),density(index_12),'ob','Linewidth',1.5)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12)
plot(feed_range,m,'--r','linewidth',1.5)
legend('Measured','Predicted mean')
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r')
alpha(0.2)
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(150,15,'Spindle speed = 3,000 (RPM)','Interpreter','Latex')
box on


subplot('Position',[0.05*3+0.8/3*2+0.02 0.2 0.8/3 0.75]);
hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_13);
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','b');
plot(feed_range,m,'--b','linewidth',1.5)
% Plot the model trained by Jinkyoo
plot(feed(index_13),density(index_13),'ob','Linewidth',1.5)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13)
plot(feed_range,m,'--r','linewidth',1.5)
legend('Measured','Predicted mean')
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r')
alpha(0.2)
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(150,15,'Spindle speed = 4,500 (RPM)','Interpreter','Latex')
box on
set(gcf,'position',[100,200,1500,300])












%first blind part
XXX = XX_measured_b{2}
YYY = YY_measured_b{2};

feed_range = 0:1:1000;
feed = XXX(:,1);
spindle = XXX(:,2);
depth_cut = XXX(:,3);
cut_direction = XXX(:,4:7);
cut_method = XXX(:,8:10);
density = YYY;

%cut == 1          input(:,4)        input(:,5)         input(:,3)    input(:,2)     input(:,2)    input(:,12)
index_11 = find(cut_direction(:,2)==1 & cut_method(:,1)==1 & depth_cut == 1 & spindle >1600 & spindle <1800 );
index_12 = find(cut_direction(:,2)==1 & cut_method(:,1)==1 & depth_cut == 1 & spindle >2700 & spindle <2900 );
index_13 = find(cut_direction(:,2)==1 & cut_method(:,1)==1 & depth_cut == 1 & spindle >4200 & spindle <4400 );

%feed            spindle                          depth_cut                     cut_direction                   cut_method
XX_11 = [feed_range',1700*ones(length(feed_range),1), 1*ones(length(feed_range),1),repmat([0 1 0 0],length(feed_range),1), repmat([1 0 0],length(feed_range),1)];
XX_12 = [feed_range',2800*ones(length(feed_range),1), 1*ones(length(feed_range),1),repmat([0 1 0 0],length(feed_range),1), repmat([1 0 0],length(feed_range),1)];
XX_13 = [feed_range',4300*ones(length(feed_range),1), 1*ones(length(feed_range),1),repmat([0 1 0 0],length(feed_range),1), repmat([1 0 0],length(feed_range),1)];

% We need to convert XX_11, XX_12, XX_13 to the five feature form used by the PMML Model
XPMML_11 = [XX_11(:,1:3), repmat([2 1],nrows,1)];
XPMML_12 = [XX_12(:,1:3), repmat([2 1],nrows,1)];
XPMML_13 = [XX_13(:,1:3), repmat([2 1],nrows,1)];


figure(2)
subplot('Position',[0.05 0.2 0.8/3 0.75]);
hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_11);
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','b');
plot(feed_range,m,'--b','linewidth',1.5)
% Plot the model trained by Jinkyoo
plot(feed(index_11),density(index_11),'ob','Linewidth',1.5)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11)
plot(feed_range,m,'--r','linewidth',1.5)
legend('Measured','Predicted mean')
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r')
alpha(0.2)
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(150,15,'Spindle speed = 1,700 (RPM)','Interpreter','Latex')
box on

subplot('Position',[0.05*2+0.8/3+0.01 0.2 0.8/3 0.75]);
hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_12);
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','b');
plot(feed_range,m,'--b','linewidth',1.5)
% Plot the model trained by Jinkyoo
plot(feed(index_12),density(index_12),'ob','Linewidth',1.5)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12)
plot(feed_range,m,'--r','linewidth',1.5)
legend('Measured','Predicted mean')
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r')
alpha(0.2)
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(150,15,'Spindle speed = 2,800 (RPM)','Interpreter','Latex')
box on


subplot('Position',[0.05*3+0.8/3*2+0.02 0.2 0.8/3 0.75]);
hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_13);
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','b');
plot(feed_range,m,'--b','linewidth',1.5)
% Plot the model trained by Jinkyoo
plot(feed(index_13),density(index_13),'ob','Linewidth',1.5)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13)
plot(feed_range,m,'--r','linewidth',1.5)
legend('Measured','Predicted mean')
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r')
alpha(0.2)
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(150,15,'Spindle speed = 4,300 (RPM)','Interpreter','Latex')
box on
set(gcf,'position',[200,200,1500,300])










%first blind part
XXX = XX_measured_b{3}
YYY = YY_measured_b{3};

feed_range = 0:1:1000;
feed = XXX(:,1);
spindle = XXX(:,2);
depth_cut = XXX(:,3);
cut_direction = XXX(:,4:7);
cut_method = XXX(:,8:10);
density = YYY;

%cut == 1          input(:,4)        input(:,5)         input(:,3)    input(:,2)     input(:,2)    input(:,12)
index_11 = find(cut_direction(:,2)==1 & cut_method(:,1)==1 & depth_cut == 1 & spindle >2030 & spindle <2230 );
index_12 = find(cut_direction(:,2)==1 & cut_method(:,1)==1 & depth_cut == 1 & spindle >2300 & spindle <2500 );
index_13 = find(cut_direction(:,2)==1 & cut_method(:,1)==1 & depth_cut == 1 & spindle >3650 & spindle <3850 );

%feed            spindle                          depth_cut                     cut_direction                   cut_method
XX_11 = [feed_range',2130*ones(length(feed_range),1), 1*ones(length(feed_range),1),repmat([0 1 0 0],length(feed_range),1), repmat([1 0 0],length(feed_range),1)];
XX_12 = [feed_range',2400*ones(length(feed_range),1), 1*ones(length(feed_range),1),repmat([0 1 0 0],length(feed_range),1), repmat([1 0 0],length(feed_range),1)];
XX_13 = [feed_range',3750*ones(length(feed_range),1), 1*ones(length(feed_range),1),repmat([0 1 0 0],length(feed_range),1), repmat([1 0 0],length(feed_range),1)];

% We need to convert XX_11, XX_12, XX_13 to the five feature form used by the PMML Model
XPMML_11 = [XX_11(:,1:3), repmat([2 1],nrows,1)];
XPMML_12 = [XX_12(:,1:3), repmat([2 1],nrows,1)];
XPMML_13 = [XX_13(:,1:3), repmat([2 1],nrows,1)];

figure(3)
subplot('Position',[0.05 0.2 0.8/3 0.75]);
hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_12);
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','b');
plot(feed_range,m,'--b','linewidth',1.5)
% Plot the model trained by Jinkyoo
plot(feed(index_11),density(index_11),'ob','Linewidth',1.5)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11)
plot(feed_range,m,'--r','linewidth',1.5)
legend('Measured','Predicted mean')
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r')
alpha(0.2)
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(150,15,'Spindle speed = 2,130 (RPM)','Interpreter','Latex')
box on

subplot('Position',[0.05*2+0.8/3+0.01 0.2 0.8/3 0.75]);
hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_12);
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','b');
plot(feed_range,m,'--b','linewidth',1.5)
% Plot the model trained by Jinkyoo
plot(feed(index_12),density(index_12),'ob','Linewidth',1.5)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12)
plot(feed_range,m,'--r','linewidth',1.5)
legend('Measured','Predicted mean')
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r')
alpha(0.2)
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(150,15,'Spindle speed = 2,400 (RPM)','Interpreter','Latex')
box on


subplot('Position',[0.05*3+0.8/3*2+0.02 0.2 0.8/3 0.75]);
hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_13);
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','b');
plot(feed_range,m,'--b','linewidth',1.5)
% Plot the model trained by Jinkyoo
plot(feed(index_13),density(index_13),'ob','Linewidth',1.5)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13)
plot(feed_range,m,'--r','linewidth',1.5)
legend('Measured','Predicted mean')
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r')
alpha(0.2)
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(150,15,'Spindle speed = 3,750','Interpreter','Latex')
box on
set(gcf,'position',[300,200,1500,300])

% Plot a single graph for the PMML report
% We choose to use a spindle speed of 2400rpm as predicted mean/variance 
% aligns well with the measured results. This is not a subplot
figure; hold on;
% Plot the model from the PMML file
[m,s] = model.score(XPMML_12);
plot(feed(index_12),density(index_12),'ob','Linewidth',1.5)
plot(feed_range,m,'--r','linewidth',1.5)
plot_variance(feed_range,(m-1.96*sqrt(s))',(m+1.96*sqrt(s))','r'); alpha(0.2);
% Plot the model trained by Jinkyoo
legend('Measured','Predicted mean')
xlabel('Feed rate $x_1$ (mm/min)','Interpreter','Latex')
ylabel('Energy density (J/mm)','Interpreter','Latex')
xlim([0,500])
ylim([0,20])
text(170,12,'Spindle speed = 2,400 (RPM)','Interpreter','Latex','fontSize',15);
set(gca,'fontSize',15)
set(gca,'XTick',0:100:500);
box on
