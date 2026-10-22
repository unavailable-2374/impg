#!/bin/bash
#SBATCH -J impg-memory-verify
#SBATCH -p skx
#SBATCH -N 1 -n 1 -c 48
#SBATCH --mem=160G
#SBATCH -t 06:00:00
#SBATCH -A TG-MCB140147
#SBATCH -o /scratch/10779/shuocao2374/tool/impg-memory/logs/production-batch.%j.log
set -euo pipefail
cd /scratch/10779/shuocao2374/tool/impg-memory
NEW=validation/bin/final
OLD=/scratch/10779/shuocao2374/tool/bin-frozen/impg-3351b8fd
EXPECTED=78244706b928055068d2bfce8072d801ba82d86206a66d66c1132f19ee3e054c
[ "$(sha256sum "$NEW" | cut -d' ' -f1)" = "$EXPECTED" ]
echo "job=$SLURM_JOB_ID host=$(hostname) started=$(date -Is)"
new_label=production-batch-new-${SLURM_JOB_ID}
old_label=production-batch-old-${SLURM_JOB_ID}
snapshot=validation/production-snapshot-${SLURM_JOB_ID}.journal
[ ! -e "$snapshot" ]
cp /scratch/10779/shuocao2374/VGP/ILS/partition/run_w5k/partition.journal "$snapshot"
chmod a-w "$snapshot"
python3 scripts/probe-production-memory.py "$NEW" "validation/$new_label" --queries 100 --threads 12 --source-journal "$snapshot" > "logs/$new_label-probe.log" 2>&1 &
new_pid=$!
python3 scripts/probe-production-memory.py "$OLD" "validation/$old_label" --queries 100 --threads 12 --source-journal "$snapshot" > "logs/$old_label-probe.log" 2>&1 &
old_pid=$!
# Each Python runner owns and cleans up its child. Failures are not passes.
status=0
wait "$new_pid" || status=1
wait "$old_pid" || status=1
[ "$status" = 0 ] || { echo 'FAIL: a production probe did not finish'; exit 1; }
python3 scripts/compare-production-probes.py "$new_label" "$old_label"
echo "PASS production equivalence and continuation; finished=$(date -Is)"
