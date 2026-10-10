# Running on a SLURM cluster (optional)

The CLI runs anywhere; on a cluster, run one job per plate. Example for the EMBL cluster layout — adapt partitions,
account and paths.

```bash
# once: environment
conda create -n diffquant python=3.11 -y && conda activate diffquant
pip install "diffquant[seg] @ git+https://github.com/<user>/diffquant"

# list plates -> one array task per plate
N=$(diffquant plates my.yaml | wc -l)
```

**Segmentation (GPU)** — `seg.sbatch`
```bash
#!/bin/bash
#SBATCH -p gpu -A <account> --gres=gpu:1 -c 4 --mem 32G -t 4:00:00
#SBATCH --array=1-N            # replace N
PLATE=$(diffquant plates my.yaml | sed -n "${SLURM_ARRAY_TASK_ID}p")
diffquant segment my.yaml --plate "$PLATE"
```

**Measurement and thresholds (CPU)** — `quant.sbatch`
```bash
#!/bin/bash
#SBATCH -p htc -A <account> -c 8 --mem 48G -t 2:00:00
#SBATCH --array=1-N
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1      # avoids BLAS thread oversubscription on many-core nodes
PLATE=$(diffquant plates my.yaml | sed -n "${SLURM_ARRAY_TASK_ID}p")
diffquant measure   my.yaml --plate "$PLATE" --jobs 8
diffquant threshold my.yaml --plate "$PLATE"
diffquant focus     my.yaml --plate "$PLATE" --jobs 8
```

```bash
J1=$(sbatch --parsable seg.sbatch)
J2=$(sbatch --parsable --dependency=afterok:$J1 quant.sbatch)
sbatch --dependency=afterok:$J2 -p htc -A <account> -c 1 --mem 8G --wrap "diffquant summarize my.yaml"
```

Notes
- Typical cost: cpsam segmentation took ~3 s per 2160 × 2160 px site on an NVIDIA L40S GPU.
- Re-running `measure` / `threshold` overwrites that plate's tables; `segment` skips sites that already have masks
  (`--overwrite` re-segments them).
- Write outputs to project storage, not to `/scratch`, if they must be kept.
