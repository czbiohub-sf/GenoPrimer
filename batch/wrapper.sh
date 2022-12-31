###################################
# Need to specify the output dir  #
# must use the name of the genome #
outdir="ensembl_GRCh38_latest"    #
###################################

echo "removing the output directory"
if [ -d $outdir ]
then 
    rm -rf $outdir
fi
mkdir $outdir

echo "removing directory slurm.out"
if [ -d "slurm.out" ]
then 
    rm -rf slurm.out
fi
mkdir slurm.out

sbatch batch_GenoPrimer.sh