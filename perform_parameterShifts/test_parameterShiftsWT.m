load ../Scer/NETseq_geneNames.mat;
load ../compare_paramsAllGenes/NETseqBestFitBinned.mat;
load ../compare_paramsAllGenes/NETseqGeneParamsBinned.mat;

model_idx = 7;

%% Make files of parameters to test
double_initFile = fopen('WT_NETseqDoubleInitParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==1
            fprintf(double_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*2);
        else
            fprintf(double_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(double_initFile,'\n');
end
fclose(double_initFile);

%%
half_initFile = fopen('WT_NETseqHalfInitParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==1
            fprintf(half_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/2);
        else
            fprintf(half_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(half_initFile,'\n');
end
fclose(half_initFile);
%%
double_elongFile = fopen('WT_NETseqDoubleElongParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==2
            fprintf(double_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*2);
        else
            fprintf(double_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(double_elongFile,'\n');
end
fclose(double_elongFile);
%%
half_elongFile = fopen('WT_NETseqHalfElongParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==2
            fprintf(half_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/2);
        else
            fprintf(half_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(half_elongFile,'\n');
end
fclose(half_elongFile);

%% Make files of parameters to test
triple_initFile = fopen('WT_NETseqTripleInitParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==1
            fprintf(triple_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*3);
        else
            fprintf(triple_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(triple_initFile,'\n');
end
fclose(triple_initFile);

%%
third_initFile = fopen('WT_NETseqThirdInitParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==1
            fprintf(third_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/3);
        else
            fprintf(third_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(third_initFile,'\n');
end
fclose(third_initFile);
%%
triple_elongFile = fopen('WT_NETseqTripleElongParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==2
            fprintf(triple_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*3);
        else
            fprintf(triple_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(triple_elongFile,'\n');
end
fclose(triple_elongFile);
%%
third_elongFile = fopen('WT_NETseqThirdElongParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==2
            fprintf(third_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/3);
        else
            fprintf(third_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(third_elongFile,'\n');
end
fclose(third_elongFile);
%% Make files of parameters to test
quadruple_initFile = fopen('WT_NETseqQuadrupleInitParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==1
            fprintf(quadruple_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*4);
        else
            fprintf(quadruple_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(quadruple_initFile,'\n');
end
fclose(quadruple_initFile);

%%
quarter_initFile = fopen('WT_NETseqQuarterInitParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==1
            fprintf(quarter_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/4);
        else
            fprintf(quarter_initFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(quarter_initFile,'\n');
end
fclose(quarter_initFile);
%%
quadruple_elongFile = fopen('WT_NETseqQuadrupleElongParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==2
            fprintf(quadruple_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*4);
        else
            fprintf(quadruple_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(quadruple_elongFile,'\n');
end
fclose(quadruple_elongFile);
%%
quarter_elongFile = fopen('WT_NETseqQuarterElongParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==2
            fprintf(quarter_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/4);
        else
            fprintf(quarter_elongFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(quarter_elongFile,'\n');
end
fclose(quarter_elongFile);

%%
no_backtrackReleaseFile = fopen('WT_NETseqNoBacktrackReleaseParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==6
            fprintf(no_backtrackReleaseFile,'%f\t',0.0);
        elseif i==10
            fprintf(no_backtrackReleaseFile,'%f\t',0.0);
        else
            fprintf(no_backtrackReleaseFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(no_backtrackReleaseFile,'\n');
end
fclose(no_backtrackReleaseFile);
%%
no_backtrackingFile = fopen('WT_NETseqNoBacktrackingParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if i==5
            fprintf(no_backtrackingFile,'%f\t',0.0);
        elseif i==9
            fprintf(no_backtrackingFile,'%f\t',0.0);
        else
            fprintf(no_backtrackingFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(no_backtrackingFile,'\n');
end
fclose(no_backtrackingFile);

