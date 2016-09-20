function [trainedClassifier, validationAccuracy] = trainClassifier(datasetTable)
% Convert input to table
datasetTable = table(datasetTable);
datasetTable.Properties.VariableNames = {'column'};
% Split matrices in the input table into vectors
datasetTable.column_1 = datasetTable.column(:,1);
datasetTable.column_2 = datasetTable.column(:,2);
datasetTable.column_3 = datasetTable.column(:,3);
datasetTable.column_4 = datasetTable.column(:,4);
datasetTable.column_5 = datasetTable.column(:,5);
datasetTable.column_6 = datasetTable.column(:,6);
datasetTable.column_7 = datasetTable.column(:,7);
datasetTable.column_8 = datasetTable.column(:,8);
datasetTable.column_9 = datasetTable.column(:,9);
datasetTable.column_10 = datasetTable.column(:,10);
datasetTable.column_11 = datasetTable.column(:,11);
datasetTable.column_12 = datasetTable.column(:,12);
datasetTable.column_13 = datasetTable.column(:,13);
datasetTable.column_14 = datasetTable.column(:,14);
datasetTable.column_15 = datasetTable.column(:,15);
datasetTable.column_16 = datasetTable.column(:,16);
datasetTable.column_17 = datasetTable.column(:,17);
datasetTable.column_18 = datasetTable.column(:,18);
datasetTable.column_19 = datasetTable.column(:,19);
datasetTable.column_20 = datasetTable.column(:,20);
datasetTable.column_21 = datasetTable.column(:,21);
datasetTable.column_22 = datasetTable.column(:,22);
datasetTable.column_23 = datasetTable.column(:,23);
datasetTable.column_24 = datasetTable.column(:,24);
datasetTable.column_25 = datasetTable.column(:,25);
datasetTable.column_26 = datasetTable.column(:,26);
datasetTable.column_27 = datasetTable.column(:,27);
datasetTable.column_28 = datasetTable.column(:,28);
datasetTable.column_29 = datasetTable.column(:,29);
datasetTable.column_30 = datasetTable.column(:,30);
datasetTable.column_31 = datasetTable.column(:,31);
datasetTable.column_32 = datasetTable.column(:,32);
datasetTable.column_33 = datasetTable.column(:,33);
datasetTable.column_34 = datasetTable.column(:,34);
datasetTable.column_35 = datasetTable.column(:,35);
datasetTable.column_36 = datasetTable.column(:,36);
datasetTable.column_37 = datasetTable.column(:,37);
datasetTable.column_38 = datasetTable.column(:,38);
datasetTable.column_39 = datasetTable.column(:,39);
datasetTable.column_40 = datasetTable.column(:,40);
datasetTable.column_41 = datasetTable.column(:,41);
datasetTable.column_42 = datasetTable.column(:,42);
datasetTable.column_43 = datasetTable.column(:,43);
datasetTable.column_44 = datasetTable.column(:,44);
datasetTable.column_45 = datasetTable.column(:,45);
datasetTable.column_46 = datasetTable.column(:,46);
datasetTable.column_47 = datasetTable.column(:,47);
datasetTable.column_48 = datasetTable.column(:,48);
datasetTable.column_49 = datasetTable.column(:,49);
datasetTable.column_50 = datasetTable.column(:,50);
datasetTable.column_51 = datasetTable.column(:,51);
datasetTable.column_52 = datasetTable.column(:,52);
datasetTable.column_53 = datasetTable.column(:,53);
datasetTable.column_54 = datasetTable.column(:,54);
datasetTable.column_55 = datasetTable.column(:,55);
datasetTable.column_56 = datasetTable.column(:,56);
datasetTable.column_57 = datasetTable.column(:,57);
datasetTable.column_58 = datasetTable.column(:,58);
datasetTable.column_59 = datasetTable.column(:,59);
datasetTable.column_60 = datasetTable.column(:,60);
datasetTable.column_61 = datasetTable.column(:,61);
datasetTable.column_62 = datasetTable.column(:,62);
datasetTable.column_63 = datasetTable.column(:,63);
datasetTable.column_64 = datasetTable.column(:,64);
datasetTable.column_65 = datasetTable.column(:,65);
datasetTable.column_66 = datasetTable.column(:,66);
datasetTable.column_67 = datasetTable.column(:,67);
datasetTable.column_68 = datasetTable.column(:,68);
datasetTable.column_69 = datasetTable.column(:,69);
datasetTable.column_70 = datasetTable.column(:,70);
datasetTable.column_71 = datasetTable.column(:,71);
datasetTable.column_72 = datasetTable.column(:,72);
datasetTable.column_73 = datasetTable.column(:,73);
datasetTable.column_74 = datasetTable.column(:,74);
datasetTable.column_75 = datasetTable.column(:,75);
datasetTable.column_76 = datasetTable.column(:,76);
datasetTable.column_77 = datasetTable.column(:,77);
datasetTable.column_78 = datasetTable.column(:,78);
datasetTable.column_79 = datasetTable.column(:,79);
datasetTable.column_80 = datasetTable.column(:,80);
datasetTable.column_81 = datasetTable.column(:,81);
datasetTable.column_82 = datasetTable.column(:,82);
datasetTable.column_83 = datasetTable.column(:,83);
datasetTable.column_84 = datasetTable.column(:,84);
datasetTable.column_85 = datasetTable.column(:,85);
datasetTable.column_86 = datasetTable.column(:,86);
datasetTable.column_87 = datasetTable.column(:,87);
datasetTable.column_88 = datasetTable.column(:,88);
datasetTable.column_89 = datasetTable.column(:,89);
datasetTable.column = [];
% Extract predictors and response
predictorNames = {'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89'};
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.column_1;
% Train a classifier
template = templateSVM('KernelFunction', 'linear', 'PolynomialOrder', [], 'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1);
trainedClassifier = fitcecoc(predictors, response, 'Learners', template, 'Coding', 'onevsone', 'PredictorNames', {'column_2' 'column_3' 'column_4' 'column_5' 'column_6' 'column_7' 'column_8' 'column_9' 'column_10' 'column_11' 'column_12' 'column_13' 'column_14' 'column_15' 'column_16' 'column_17' 'column_18' 'column_19' 'column_20' 'column_21' 'column_22' 'column_23' 'column_24' 'column_25' 'column_26' 'column_27' 'column_28' 'column_29' 'column_30' 'column_31' 'column_32' 'column_33' 'column_34' 'column_35' 'column_36' 'column_37' 'column_38' 'column_39' 'column_40' 'column_41' 'column_42' 'column_43' 'column_44' 'column_45' 'column_46' 'column_47' 'column_48' 'column_49' 'column_50' 'column_51' 'column_52' 'column_53' 'column_54' 'column_55' 'column_56' 'column_57' 'column_58' 'column_59' 'column_60' 'column_61' 'column_62' 'column_63' 'column_64' 'column_65' 'column_66' 'column_67' 'column_68' 'column_69' 'column_70' 'column_71' 'column_72' 'column_73' 'column_74' 'column_75' 'column_76' 'column_77' 'column_78' 'column_79' 'column_80' 'column_81' 'column_82' 'column_83' 'column_84' 'column_85' 'column_86' 'column_87' 'column_88' 'column_89'}, 'ResponseName', 'column_1', 'ClassNames', [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14]);

% Perform cross-validation
partitionedModel = crossval(trainedClassifier, 'KFold', 5);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

%% Uncomment this section to compute validation predictions and scores:
% % Compute validation predictions and scores
% [validationPredictions, validationScores] = kfoldPredict(partitionedModel);