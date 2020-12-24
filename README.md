## Comprehensive analysis and GWAS of biomass, chlorophyll, seed and salinity tolerance related traits in rice :ear_of_rice: :paw_prints:
### Section 1: phenotype processing
- Notebooks/Phenotypes_compilation.ipynb :arrow_right: Prepares working genotype means from plant data
- Notebooks/Phenotype_stats.ipynb :arrow_right: Basic analytics on phenotypes i.e. histograms, density plots, Shapiro-tests
- Notebooks/Broad-sense_heritability.ipynb :arrow_right: Calculates broad-sense heritability from total and genotype variance
- R scripts/Random_effects_modelling_for_heritability.R :arrow_right: Estimates trait heritability by modelling genotype, condition and their interactions as random effects
- R scripts/Marker-based_heritability.R :arrow_right: Estimates heritability from genomic kinship

### Section 2: genotype preprocessing
- Shell scripts/1.download_176vcf_data.sh :arrow_right: Uses curl to download individual VCF files from the 3000-rice genome project
- Shell scripts/2.gzip_to_bgzip.sh :arrow_right: Converts gziped VCF filed to bgzip compression
- Shell scripts/3.combine_vcf.sh :arrow_right: Combines individual VCF files into one
- Shell scripts/4.beagle_imputation.sh :arrow_right: Imputes missing marker genotypes using Beagle 5.1
- Notebooks/Imputation_accuracy.ipynb :arrow_right: Assessment of imputation accuracy
- Python scripts/Make_working_files.py :arrow_right: Prepares a number of working files
- Python scripts/Make_hmp.py :arrow_right: Prepares a hapmap genotype file
- Shell scripts/5.plink_conversion_and_pruning.sh :arrow_right: Prepares plink files and estimates effective number of markers

### Section 3: genomic predictions, phenotype transformations
- R scripts/Genomic_predictions.R :arrow_right: Uses ridge regression in mixed.solve() to predict phenotypes
- Python scripts/Transformations_p1.py :arrow_right: Prepares a shell script for WarpedLMM transformation
- Shell scripts/7.transform_phenotypes.sh :arrow_right: Executes WarpedLMM
- Python scripts/Transformations_p2.py :arrow_right: Compiles WarpedLMM results

### Section 4: population structure and GWAS
- R scripts/Population_structure_estimation.R :arrow_right: Population structure estimation using genomic scatter plots, PCA and k-means clustering
- Shell scripts/6.fastStructure1-15.sh :arrow_right: Employs fastStructure for population structure estimation and finds appropriate number of subpopulations
- Python scripts/Split_populations.py :arrow_right: Splits working files into subpopulations according to population structure
- R scripts/GAPIT_for_GWAS.R :arrow_right: Tests markers for phenotype association using the BLINK algorithm and CMLM
- R scripts/9.LD_decay.sh :arrow_right: Determines extent of linkage disequilibrium

### Section 5: downstream/post-GWAS analytics
- R scripts/Plotting_GWAS_results.R :arrow_right: Prepares manhattan and quantile-quantile plots
- Shell scripts/8.blast.sh :arrow_right: BLAST for finding physical locations and ranges of known genes
- Notebooks/Significant_Intergenic_markers.ipynb :arrow_right: Compiles significant and suggestive marker associations from GWAS that are within known gene regions
- R scripts/Dendogram_and_second_gene_expression_heatmap.R :arrow_right: Clusters genes by dendograms and heatmaps
- Notebooks/Slice_VCF.ipynb :arrow_right: Extracts intergenic markers
- R scripts/Beautiful_Exon_Extractor.R :arrow_right: Extracts exons from pairs of genes and CDSes
- R scripts/Beautiful_Intron_Masker.R :arrow_right: Masks introns from gene-CDS pairs
- Notebooks/SNP_effects_and_haplotype_testing.ipynb :arrow_right: Deciphers protein level consequences of polymorphisms and tests alleles by ANOVA and Student's t test
- Notebooks/Multiple_testing_correction_and_LD_statistics.ipynb :arrow_right: Calculates FDR-adjusted p values using the Benjamini-Hochberg method. Evaluates LD for markers and QTLs
- Notebooks/Plots.ipynb :arrow_right: Miscellaneous visualizations

