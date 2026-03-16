#!/bin/bash
#SBATCH --job-name=phase8a_long
#SBATCH --partition=compute2
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=2-12:00:00
#SBATCH --output=/work/vmo703/tristan-mp-v2/logs/%x_%j/slurm.out
#SBATCH --error=/work/vmo703/tristan-mp-v2/logs/%x_%j/slurm.err

set -euo pipefail

source /apps/hpcx/2.18.1/hpcx-init.sh
hpcx_load

LOGDIR="/work/vmo703/tristan-mp-v2/logs/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir -p "$LOGDIR"

INPUT="/work/vmo703/tristan-mp-v2/runs/phase8/input_phase8a_long_control"
OUTPUT="/work/vmo703/tristan-mp-v2/runs/phase8/out_phase8a_long_control"
BIN="/work/vmo703/tristan-mp-v2/build_3d_8a/src/tristan-v2.xc"

echo "=============================================="
echo " phase 8A-long source-plus-buffer control"
echo " log directory:"
echo "   $LOGDIR"
echo " input:  $INPUT"
echo " output: $OUTPUT"
echo " binary: $BIN"
echo "=============================================="
echo ""

[[ -f "$INPUT" ]] || { echo "ERROR: missing input file: $INPUT"; exit 2; }

mkdir -p "$OUTPUT"

export OPAL_PREFIX=/apps/hpcx/2.18.1/ompi
export LD_LIBRARY_PATH=/work/vmo703/hdf5-parallel/lib:/apps/anaconda3/2024.02-1/lib:$LD_LIBRARY_PATH
export HDF5_USE_FILE_LOCKING=FALSE

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp
export OMPI_MCA_coll_ucc_enable=0
export OMPI_MCA_coll_hcoll_enable=0
export UCX_TLS=sm,self,tcp
export UCX_NET_DEVICES=lo
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

MPI=mpirun

$MPI --mca pml ob1 --mca btl self,tcp \
     --mca rmaps_base_oversubscribe 1 \
     -np $SLURM_NTASKS \
     $BIN -i "$INPUT" -o "$OUTPUT"

echo ""
echo "phase 8A-long run completed."
echo "logs available at: $LOGDIR"
