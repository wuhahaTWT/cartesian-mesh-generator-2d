#!/bin/bash
# Mature solver only; bash resource monitoring avoids adding Python to its image.
case_directory=$1
budget=$2
cd "$case_directory" || exit 2
source /usr/lib/openfoam/openfoam2606/etc/bashrc >/dev/null 2>&1
checkMesh > checkMesh.log 2>&1
checkMesh -allGeometry -allTopology > checkMesh-extended.log 2>&1
TIMEFORMAT='{"wallSeconds":%0R,"userSeconds":%0U,"systemSeconds":%0S}'
{ time timeout --signal=TERM --kill-after=10s "$budget" simpleFoam > solve.log 2>&1; } 2> inner-time.json &
solve_pid=$!
maximum=0
while kill -0 "$solve_pid" 2>/dev/null; do
    total=0
    for status in /proc/[0-9]*/status; do
        [ -r "$status" ] || continue
        while read -r key value rest; do
            if [ "$key" = 'VmRSS:' ]; then total=$((total+value)); break; fi
        done < "$status" 2>/dev/null || true
    done
    if [ "$total" -gt "$maximum" ]; then maximum=$total; fi
    printf '%s\n' "$maximum" > peak-container-process-rss-kib.txt
    sleep .5
done
wait "$solve_pid"
code=$?
printf '{"exitCode":%s,"peakProcessRssMiB":%.3f}\n' "$code" "$(printf '%s' "$maximum" | awk '{print $1/1024}')" > inner-status.json
exit "$code"
