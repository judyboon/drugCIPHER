function [Drug_Gene_Score_Matrix] = drugCIPHER_SingleS_Overall(Drug_Target_Relation, PPI_Adjacent_Matrix, Drug_Similarity_Matrix,Scale)

%%%%%%%%%%

%% Drug_Target_Relation:
%%% Drug_Targe_Relation contains targets information for each drug, every
%%% line represents a drug, and the corresponding known targets are
%%% seperated by '\t'. If the line is '-1', there is no known target for
%%% that drug. The targets are represented by the gene index in the PPI
%%% network. The index are begin from 0. 


%% PPI_Adjacent_Matrix
%%% PPI_Adjacent_Matrix is the adjacent matrix for the input PPI network.
%%% 1 presents having edge; otherwise, is 0

%% Drug_Similarity_Matrix1, Drug_Similarity_Matrix2
%%% The drug similarity matrix, corresponding to the index in Drug_Target_Relation.

%% Repeat_Times
%%% The random repeat times of the validation procedure, the default value
%%% is 100;



%% Set the default value



%% Begin Alogrithm
%%%%% 
%%%%%
%%%%%

%% Compute Shortest Distance in the PPI Network
disp('Computing The Shortest Distance Between Nodes...');


Suffix_Exame = regexp(PPI_Adjacent_Matrix,'\w*txt','match');
tag = 0;
if isempty(Suffix_Exame)
    tag = 1;
    Suffix_Exame = regexp(PPI_Adjacent_Matrix,'\w*mat','match');
    if isempty(Suffix_Exame)
        disp('error, please check the input format');
        return;
    end
end

if tag == 0
    PPI_Adjacent_Matrix_data = load(PPI_Adjacent_Matrix);
else
    Struct_Data = load(PPI_Adjacent_Matrix);
    Field_Name = fieldnames(Struct_Data);
    PPI_Adjacent_Matrix_data = getfield(Struct_Data,Field_Name{1});
end
%PPI_Adjacent_Matrix_data(find(PPI_Adjacent_Matrix_data == -1)) = 0;
Sparse_Matrix = sparse(PPI_Adjacent_Matrix_data);
clear PPI_Adjacent_Matrix_data;
PPI_Shortest_Dist = graphallshortestpaths(Sparse_Matrix,'Directed','false');
[GeneNum,t] = size(PPI_Shortest_Dist);

%% Compute Gene2drug Closeness
disp('Computing Gene to Drug Closeness Matrix...');

Suffix_Exame = regexp(Drug_Target_Relation,'\w*txt','match');
if isempty(Suffix_Exame)
    disp('error, please check the input format');
    return;
end

fid = fopen(Drug_Target_Relation);
Drug_Target_Relation_Number = 0;
DrugNum = 1;
line = fgetl(fid);
while ischar(line)
    tArray = regexp(line,'\t','split');
    if ~isempty(line)
        Array = [];
        [t,TargetNum] = size(tArray);
        if strcmp(tArray(1),'-1')
            Drug2Targets{DrugNum} = [];
        else
            Drug_Target_Relation_Number = Drug_Target_Relation_Number + TargetNum;
            Array(TargetNum) = 0;
            for i = 1:TargetNum
                Array(i) = str2num(tArray{i});
            end           
            Drug2Targets{DrugNum} = Array;            
        end
        line = fgetl(fid);
        DrugNum = DrugNum + 1;
    else
        break;
    end
end
DrugNum = DrugNum - 1;
fclose(fid);

Gene2Drug_Closeness(GeneNum,DrugNum) = 0;

for i = 1:GeneNum
    for j = 1:DrugNum
        Array = Drug2Targets{j};
        if ~isempty(Array)
            Gene2Drug_Closeness(i,j) = sum(exp(-(PPI_Shortest_Dist(i,Array+1)).^2));
        else
            Gene2Drug_Closeness(i,j) = 0;
        end
    end
end


%% Compute Drug Gene Score
disp('Loading Drug Similarity Matrixes...');

Suffix_Exame = regexp(Drug_Similarity_Matrix,'\w*txt','match');
tag = 0;
if isempty(Suffix_Exame)
    tag = 1;
    Suffix_Exame = regexp(Drug_Similarity_Matrix,'\w*mat','match');
    if isempty(Suffix_Exame)
        disp('error, please check the input format');
        return;
    end
end

if tag == 0
    Drug_Similarity_Matrix_data = load(Drug_Similarity_Matrix);
else
    Struct_Data = load(Drug_Similarity_Matrix);
    Field_Name = fieldnames(Struct_Data);
    Drug_Similarity_Matrix_data = getfield(Struct_Data,Field_Name{1});
end

Drug_Similarity_Matrix_data = Drug_Similarity_Matrix_data.^Scale;

Drug_Gene_Score_Matrix(DrugNum,GeneNum) = 0;

for i = 1:DrugNum %For all drugs, TargetsSet is the validation set; the number i-th drug
    fprintf('    Computing Target Profile for the %dth drug...\n',i);
    Drug_i_similarity = Drug_Similarity_Matrix_data(i,:)';

    Drug_Gene_Score_Matrix(i,:) = corr(Drug_i_similarity,Gene2Drug_Closeness');
    
    %%%%%%%%%%%%%%%%
end


