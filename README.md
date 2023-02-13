# SARS2EVO

The site includes all custom scripts and files used in the study of “**Underlying driving forces of the SARS-CoV-2 evolution: immune evasion and ACE2 binding affinity**”.

The code includes python scripts, C++ scripts, and C# scripts. Both Linux and Windows platforms are required. We recommend running the scripts on a Windows working station with 128Gb RAM and enabling the Windows subsystem Linux (WSL) system to avoid intermediate file transfer which can be very time-consuming. In this documentation, unless otherwise noted, all steps are run under the Windows system.

Most of the scripts are written in Microsoft C# language based on the .net framework, the original code and project files are provided. To run the scripts, [Microsoft Visual Studio](https://visualstudio.microsoft.com/) is needed. The directory of the input and output should be named in the script before hitting the Compile and run button.

# Data collection

A total number of 6,484,070 high-quality open-access SARS-CoV-2 sequences and corresponding metadata were downloaded from the [UShER website](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/) on 23 November 2022 as the input. The following files were downloaded and unzipped:

_public-latest.all.masked.pb.gz_

_public-latest.metadata.tsv.gz_

_public-latest.all.masked.vcf.gz_

(http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/)

# Calculate mutation rate (Figure 1)

**2.1** **Obtain mutation and metadata for each sequence**

2.1.1 Transfer the unzipped VCF file to _mut5col_ format using **Vcf2mut5col** (Linux).

*“Vcf2mut5col.out _public-latest.all.masked.vcf_  [output directory]”**

*This will generate a _.mut5col_ file.

2.1.2 Transfer the _.mut5col_ file to _sample_mutlist.tsv_ using **Mutresult2Sample_mutlistLINUX** (Linux).

*“Mut5col2Sample_mutlistLINUX.out [directory]”**

*Both input and output will be in this directory. The _sample_mutlist.tsv_ is a txt file that stores mutation information of each sequence in the format of [SequenceName]\tMutation1 Mutation2 ... MutationN.

**2.2** **Sequence sampling**

Random selection of 200 sequences each month using **SampleMutlistFilterAndRandomPick**. The input is the _sample_mutlist.tsv_ file and the  _public-latest.metadata.tsv_  file_.This will generate the _sampled_sample_mutlist.tsv_ and the _sampled_metadata.tsv_  file.

**2.3** **Sequence mutation annotation**

Before plotting, the _sampled_sample_mutlist.tsv file_ was further annotated, transferring each nucleotide mutation to amino acid mutation, generating a result table file _sampled_sample_mutlist.transfer.tsv_ using **SampleMutlistTransfer**.

*The _sampled_sample_mutlist.transfer.tsv_ can be directly combined with the  _sampled_metadata.tsv_  file in the [Excel](https://www.microsoft.com/zh-cn/microsoft-365/?omkt=en-US&rtc=1) for they share the same sequence order. 

**2.4** **Mutation rate calculate**

After combining, the transferred table was then plotted using the python script **segments_fit.py**  (for the automatic piecewise linear regression ) and the **[R ggplot](https://www.rdocumentation.org/packages/ggplot2/versions/3.4.1)** function (for simple linear regression).

# Identification of all mutation events

The downloaded [protobuf](https://protobuf.dev/) file _public-latest.all.masked.pb.gz_  was first transferred to a readable JSON file using [matUtils](https://usher-wiki.readthedocs.io/en/latest/QuickStart.html#matutils) (Linux) from the UshER toolkit. (https://usher-wiki.readthedocs.io/en/latest/)

_“matUtils extract -i ./public-latest.all.masked.pb.gz -d [_YOUR DIRECTORY_] -j [_MAT1128_].json”_

*This step should be run on Linux and can require memory usage bigger than 65Gb.

Then we search for all mutation events from the MAT tree using **Json2MutationEvent_C#**. It will generate a _.json.mutevent_ file. The file contains all filtered sequences (see methods section) and their additional mutations relative to ancestors.

*The input files are the _public-latest.metadata.tsv_ and the _.json_ file. 

# Calculate the mutation distribution on the genome (Figure 1)

Mutation distribution on the genome and collection date was extracted from the obtained _.mutevent_ file using **GenomeMutationDistribution**.

# Calculate the mutation incidence in different SARS-CoV-2 lineages (Figure 2)

The mutation incidence of the most frequent mutations in different SARS-CoV-2 lineages was calculated using **FindCladeSpecificMutation** which can automatically generate several heatmap tables for each gene.

*The input file is the *.json.mutevent* and the output files are heatmap.txt files:

The heatmap was then plotted with the **R pheatmap** package. The principal component analysis (PCA) plot was plotted with the **[R fviz](https://www.rdocumentation.org/packages/factoextra/versions/1.0.7/topics/fviz_pca)** package.

# Calculate the Silhouette Coefficient (Figure 2)

After obtaining the top most frequent mutations table, we calculated the Silhouette Coefficient of each gene and clade using **CalSilhouetteCoefficient**.

*The input files are the _CladeHeatmap.txt_ file of each gene.

# Calculate the distribution of escape mutations (Figure 3)

**7.1** **Calculate the escape score**

The antibody spectrum, neutralizing activity, antibody epitope group, and raw mutation escape score were obtained from previous studies (see methods). We combined the raw mutation escape score with antibody-neutralizing activity using **CalMutEscapeScore**.

*The input DMS and output files are the following:

**Files**
use_abs_res.csv
NeutralWTBA125_Cross.txt
single_mut_effects.csv
EscapeScore_PKU.single.txt
EscapeScore_PKU.12.txt

**Description**
Raw escape scores.
Antibody neutralization data.
ACE2 binding and RBD expression data.
The processed escape score.
The processed escape score of the 12 epitopes.

**7.2** **Assign escape mutations**

The mutation that significantly reduced the affinity to any of the 12 antibody types was defined as an immune escape mutation. The escape mutations were assigned based on the escape score in the file _EscapeScore_PKU.12.txt_  using **EachCladeEscapeScoreTo12Epitope**. The output escape mutation distribution heatmap file _CladeRBDMut12Group.txt_  was also calculated using this script.

**7.3** **Calculate the antibody pressure**

The immune pressure exerted on a particular epitope region is calculated by summing the neutralizing activities of all antibodies that belong to this epitope group using **Cal12EpiPressureDistribution**.

*The input file is the antibody neutralization table _NeutralWTBA125_Cross.txt_ and the output files are the _AntibodySpectrum12_Count.txt_ and the _AntibodySpectrum12_logIC50.txt_ which represent the antibody number and pressure act on the 12 epitopes.

**7.4** **The gene set enrichment analysis (GSEA)**

The GSEA analysis was also performed using the script **FindCladeSpecificMutation**, as one of its functions.

*It will generate 2 files as output:

**File**
Mutation_ES_Curve.txt
Mutation_ES_Pvalue.txt

**Description**
The GSEA plot curve data.
The p-value after 50000 randomizations repeating.

# Calculate the evolution trajectory (Figure 4)

**8.1** **Construct a tree of Lineages**

To calculate the evolution trajectory, it is necessary to restore the evolutionary relationship of all Lineages. We obtained the Lineage list from the [Pango](https://pangolin.cog-uk.io/) website and refined the relationships into _lineageRelations.tsv_ using **LineageRelations**.
(https://github.com/cov-lineages/lineages-website)

*input: _[lineages.yml](https://github.com/cov-lineages/lineages-website/blob/master/data/lineages.yml)_

Output: _lineageRelations.tsv_

**8.2** **Summarize genome mutations and metadata for each Lineage**

We first combined the _sample_mutlist.tsv_ file with the _metadata.tsv_ file using **SeqFindMeta** (Linux). Then, we counted the sampling time and genome mutation of each Lineage sequence using **LineageMutationMetadata**, generating the _Lineage.Date.AAmut_ file.

**8.3** **Calculate the evolution trajectory**

Finally, we combined the above result with the escape score and ACE2 binding data, using **LineageTrajectory**, generating a table _LineagePlotData.txt_  that can be directly used for **[R ggplot](https://www.rdocumentation.org/packages/ggplot2/versions/3.4.1)**.

# Multivariate linear regression (Figure 4)

**9.1** **Calculation of the RoHo**

The [RoHo](https://usher-wiki.readthedocs.io/en/latest/tutorials.html#calculating-by-mutation-roho-with-matutils-summary-and-python) value of each mutation event was calculated using the [matUtils](https://usher-wiki.readthedocs.io/en/latest/matUtils.html) (Linux) from the UshER toolkit.

_“matUtils summary -i ./public-latest.all.masked.pb.gz -E RoHo.tsv”_

And then extract the corresponding data of the target clade into _MutationRoHo.txt_  using **CalMutRoHo**.

**9.2** **Multivariate linear regression**

The mutation incidence, RoHo, and corresponding escape score and ACE2 binding data were combined manually using [Microsoft Office Excel](https://www.microsoft.com/en-us/microsoft-365/excel). The multivariate linear regression was conducted by the **[R lm](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/lm)** function.

_“lm.sol = lm(EventFraction~ scale(DMSACE) + scale(DMSESC), data = RoHo_WT)”_
_“summary(lm.sol)”_
