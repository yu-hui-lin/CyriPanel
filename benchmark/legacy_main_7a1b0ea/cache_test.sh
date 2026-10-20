#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J cache_test
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -t 00:20:00
#SBATCH -o cache_test.out
#SBATCH -e cache_test.err
#SBATCH --mail-user=yhlin.md05@nycu.edu.tw
#SBATCH --mail-type=BEGIN,END,FAIL
# Does SLURM count file cache as job memory on these nodes?
# The job writes an 8 GiB file of zeros, which needs almost no process memory, waits, then deletes it.
# If file cache is counted, MaxRSS will be close to 8 GiB. If not, it stays near zero.
W=/work/u7715055/staging/biology/u7715055
T=$W/cache_test_$SLURM_JOB_ID.tmp
echo "node $(hostname)  cgroup filesystem $(stat -fc %T /sys/fs/cgroup)"
df -T "$W" | tail -1
dd if=/dev/zero of="$T" bs=1M count=8192 conv=fsync status=none
echo "written $(du -h "$T" | cut -f1) at $(date +%T)"
sleep 90
rm -f "$T"
