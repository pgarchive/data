#!/bin/bash
#SBATCH --array=1-16290
#SBATCH --job-name=CAAS_test
#SBATCH --dependency=afterany:33473911
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=06:00:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+
module load java
module load EMBOSS

#Define PATH argument
REF_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED/
RESEQUENCED_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Trees/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/FastCAAS/fastcaas/fastcaas
PROTEIN_ALIGNMENTS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_PROTEIN_ALL/Protein_alignments_filtered/
CONFIG_FILES=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/GenePhenMap/CAAS/FINAL_CAAS/Configs/Qualitative_ALLBG/subset8/
FAMILIES=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/FastCAAS/fastcaas/_tests/sp2fam.210727.tab
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/GenePhenMap/CAAS/FINAL_CAAS/Results/Phylo/Qualitative_ALLBG/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_fixed_annotations_ALL/Prepare_ref_DNA_alignment/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}SortCDS_input_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
firstchar=${gene_name:0:1}
secondchar=${gene_name:0:2}
for filepath in $(ls ${CONFIG_FILES}*.txt);
do for filepath2 in $(ls ${PROTEIN_ALIGNMENTS}/${firstchar}/${secondchar}/${gene_name}*.fasta); 
do species_name=$(echo $filepath2 | rev | cut -d'/' -f1 | rev | cut -d \. -f 2)
trait_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d '.' -f 1)
${BIN} run -t $filepath \
-a ${PROTEIN_ALIGNMENTS}/${firstchar}/${secondchar}/${gene_name}.${species_name}.filter2.phy --fmt phylip-relaxed \
--species ${FAMILIES} \
--cycles=1000 \
--scenarios=1,2 \
-j ${OUTDIR}${trait_name}/${firstchar}/${secondchar}/ \
--maxpval=1
done
done

