#!/usr/bin/env bash
# Report a Slurm job's state to Snakemake as one of: success | failed | running.
# Snakemake polls this with the job id returned by `sbatch --parsable`.
set -euo pipefail

jobid="$1"

# sacct can lag briefly after submission; fall back to scontrol, and treat an
# unknown job as still running rather than failing the whole workflow.
state="$(sacct -j "$jobid" --format=State --noheader --parsable2 2>/dev/null \
         | head -n1 | cut -d' ' -f1 || true)"

if [[ -z "$state" ]]; then
  state="$(scontrol show job "$jobid" 2>/dev/null \
           | grep -oP 'JobState=\K\S+' || true)"
fi

case "$state" in
  COMPLETED)                                  echo "success" ;;
  BOOT_FAIL|CANCELLED|DEADLINE|FAILED|NODE_FAIL|OUT_OF_MEMORY|PREEMPTED|REVOKED|TIMEOUT)
                                              echo "failed"  ;;
  *)                                          echo "running" ;;
esac
