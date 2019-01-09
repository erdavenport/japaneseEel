# 2019-01-09

# Create population map file with four metapopulations: 

# Make output directories:
mkdir data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/four_meta_populations/
mkdir data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/four_meta_populations/
mkdir data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/four_meta_populations/

# Stack depth = 3:  
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/ \
	-b 1 \
	-O data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/four_meta_populations/ \
	-M data/sample_info/population_file_for_stacks_no_JJ-107_four_meta_populations.txt \
	-W data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/whitelist_loci_after_coverage_filtering.txt \
	-t 4 \
	-s \
	-p 1 \
	-r 0.6 \
	-m 3 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--phylip \
	--verbose
	
# Stack depth = 6:  
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/ \
	-b 3 \
	-O data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/four_meta_populations/ \
	-M data/sample_info/population_file_for_stacks_no_JJ-107_four_meta_populations.txt \
	-W data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/whitelist_loci_after_coverage_filtering.txt \
	-t 4 \
	-s \
	-p 1 \
	-r 0.6 \
	-m 6 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--phylip \
	--verbose

# Stack depth = 10:  
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/ \
	-b 4 \
	-O data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/four_meta_populations/ \
	-M data/sample_info/population_file_for_stacks_no_JJ-107_four_meta_populations.txt \
	-W data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/whitelist_loci_after_coverage_filtering.txt \
	-t 4 \
	-s \
	-p 1 \
	-r 0.6 \
	-m 10 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--phylip \
	--verbose

# Recreate population map with shapes fixed:
scripts/make_map_for_eelseq_project.R

# Remake pngs for readme
sips -s format png results/1_general_info/map_sample_locations.pdf --out results/1_general_info/map_sample_locations.png
sips -s format png results/1_general_info/timeline_sample_collection.pdf --out results/1_general_info/timeline_sample_collection.png

