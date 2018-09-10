load ../Scer/spt4_NETseq_geneNames.mat;
load ../compare_paramsAllGenes/spt4_NETseqBestFitBinned.mat;
load ../compare_paramsAllGenes/spt4_NETseqGeneParamsBinned.mat;

model_idx = 7;

%% Make files of parameters to test
double_initFile = fopen('spt4_NETseqDoubleInitParams.txt','w');
for g=1:length(spt4_NETseq_geneNames)
    for i=1:size(spt4_NETseqGeneParamsBinned{model_idx},2)-3
        if i==1
            fprintf(double_initFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i)*2);
        else
            fprintf(double_initFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(double_initFile,'\n');
end
fclose(double_initFile);

%%
half_initFile = fopen('spt4_NETseqHalfInitParams.txt','w');
for g=1:length(spt4_NETseq_geneNames)
    for i=1:size(spt4_NETseqGeneParamsBinned{model_idx},2)-3
        if i==1
            fprintf(half_initFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i)/2);
        else
            fprintf(half_initFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(half_initFile,'\n');
end
fclose(half_initFile);
%%
double_elongFile = fopen('spt4_NETseqDoubleElongParams.txt','w');
for g=1:length(spt4_NETseq_geneNames)
    for i=1:size(spt4_NETseqGeneParamsBinned{model_idx},2)-3
        if i==2
            fprintf(double_elongFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i)*2);
        else
            fprintf(double_elongFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(double_elongFile,'\n');
end
fclose(double_elongFile);
%%
half_elongFile = fopen('spt4_NETseqHalfElongParams.txt','w');
for g=1:length(spt4_NETseq_geneNames)
    for i=1:size(spt4_NETseqGeneParamsBinned{model_idx},2)-3
        if i==2
            fprintf(half_elongFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i)/2);
        else
            fprintf(half_elongFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(half_elongFile,'\n');
end
fclose(half_elongFile);
%%
no_backtrackReleaseFile = fopen('spt4_NETseqNoBacktrackReleaseParams.txt','w');
for g=1:length(spt4_NETseq_geneNames)
    for i=1:size(spt4_NETseqGeneParamsBinned{model_idx},2)-3
        if i==6
            fprintf(no_backtrackReleaseFile,'%f\t',0.0);
        elseif i==10
            fprintf(no_backtrackReleaseFile,'%f\t',0.0);
        else
            fprintf(no_backtrackReleaseFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(no_backtrackReleaseFile,'\n');
end
fclose(no_backtrackReleaseFile);
%%
no_backtrackingFile = fopen('spt4_NETseqNoBacktrackingParams.txt','w');
for g=1:length(spt4_NETseq_geneNames)
    for i=1:size(spt4_NETseqGeneParamsBinned{model_idx},2)-3
        if i==5
            fprintf(no_backtrackingFile,'%f\t',0.0);
        elseif i==9
            fprintf(no_backtrackingFile,'%f\t',0.0);
        else
            fprintf(no_backtrackingFile,'%f\t',spt4_NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(no_backtrackingFile,'\n');
end
fclose(no_backtrackingFile);

