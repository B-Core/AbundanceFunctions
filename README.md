# Notes about functions in this repository

See 
 * AbundanceExampleWorkflow.Rmd 
 * analysisR.md

For examples of RNA-Seq workflows that use these functions.

## Summarize.R (updated 6 Dec., 2016) 
Summarize\_by\_some\_custom\_ID in Summarize.R of this repo. still needs some improvement because it relies on the Oligo package from BioConductor, which is not well-coordinated with Affymetrix releases. It's likely that under the hood, oligo::rma() relies on properly processing the .cdf file, which we know for a fact does not always happen with certain array types.

One perhaps insightful observation, for instance, is that summarized_mat generated from:

```
source("~/Git_repos/AbundanceFunctions/ExtractTransformLoad.R")
path2CELfiles2 = "/Users/fishema/Desktop/CEL_files_for_testing/HTA/"
microarrayProcessing_ls = processCELfiles(path2CELfiles2)
summarized_mat  = oligo::rma(microarrayProcessing_ls$MicroarrayObject, background=F, normalize=F)
```

for the human HTA array does not produce the same number of rows in its summarized matrix as summarized_mat from 

```
probeset_df = getProbeInfo(microarrayProcessing_ls$MicroarrayObject,field=c('fid','fsetid'),probeType = "pm",target = "probeset")
source("~/Git_repos/AbundanceFunctions/Summarize.R")
feature_ID_vec = probeset_df$fid
custom_ID_vecInst = probeset_df$fsetid
NormWithRownamesInst_mat = NormedMats_ls[[1]]
oligoFeatureSetObjInst = microarrayProcessing_ls$MicroarrayObject
summarized_mat = Summarize_by_some_custom_ID(oligoFeatureSetObj = oligoFeatureSetObjInst, NormWithRownames_mat = NormWithRownamesInst_mat, featureID_v = feature_ID_vec,customID_v = custom_ID_vecInst)
```

Another way of thinking about this is that probeset_df is not generating appropriate feature and probeset IDs.

We suspect/hope that this issue might be resolved once we get a tool to process/parse .cdf files.

In the meantime, Summarize\_by\_some\_custom\_ID throws a warning if the user provides oligoFeatureSetObj, NormWithRownames\_mat, featureID\_v, **and** customID\_v and the two summarized_mats don't have the same nrows().
