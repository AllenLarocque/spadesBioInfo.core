#!/usr/bin/env bash
# memory-sampler.sh <runId> <watchPid> <outFile> <interval>
# Appends one row per interval while <watchPid> is alive: PSS (anon, file) summed over
# <watchPid> and every process whose environment holds SPADES_RUN_ID=<runId>, plus the
# cgroup's memory.current. The watched PID is always included: Sys.setenv() in R does not
# show in its own /proc/<pid>/environ. Exits by itself when the watched process is gone.
runId=$1; watchPid=$2; out=$3; interval=$4

echo $$ > "$out.pid"
echo "timestamp,pssAnonKiB,pssFileKiB,nProcs,memCurrentKiB" > "$out"

cg=$(sed -n 's/^0:://p' /proc/self/cgroup)
memFile="/sys/fs/cgroup${cg}/memory.current"

while kill -0 "$watchPid" 2>/dev/null; do
  ts=$(date +%s.%N)
  pids=$({ echo "$watchPid"
           grep -lzx "SPADES_RUN_ID=$runId" /proc/[0-9]*/environ 2>/dev/null |
             sed 's#/proc/\([0-9]*\)/environ#\1#'; } | sort -un)
  files=()
  for p in $pids; do files+=("/proc/$p/smaps_rollup"); done
  # cat skips a process that vanished; nProcs counts the rollups actually read.
  sums=$(cat "${files[@]}" 2>/dev/null |
    awk '/^Pss_Anon:/ {a += $2; n++} /^Pss_File:/ {f += $2} END {printf "%d,%d,%d", a, f, n}')
  cur=""
  [ -r "$memFile" ] && cur=$(( $(cat "$memFile" 2>/dev/null || echo 0) / 1024 ))
  IFS=, read -r a f n <<< "$sums"
  echo "$ts,$a,$f,$n,$cur" >> "$out"
  sleep "$interval"
done
