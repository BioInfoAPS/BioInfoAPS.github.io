#!/bin/bash
while read line
do
        [[ `grep "BlastOutput_db" <<< $line` ]] && echo -n ">" <<< $line
        [[ `grep -E "BlastOutput_db|Hsp_hseq" <<< $line` ]] && sed -n "s:.*>\(.*\)</.*:\1:p" <<< $line
done