#PBS -lnodes=1:cores16:ppn=16
#PBS -lwalltime=12:00:00

module load mpicopy
module load openmpi

cd $HOME/EMRI/code/sims/JOB_ID/run
sed -i -- "s@RUN_DIRECTORY@$TMPDIR/JOB_ID/run@g" EMRI.param
sed -i -- "s@OUT_DIRECTORY@$TMPDIR/JOB_ID/out@g" EMRI.param

mpicopy $HOME/EMRI/code/sims/JOB_ID -o "$TMPDIR"

module load openmpi/gnu

cd "$TMPDIR"/JOB_ID/run

time mpiexec Gadget2 EMRI.param > output

tail -n1000 output > output_tail
head -n1000 output > output_head

tail -n1000 "$TMPDIR"/JOB_ID/out/cpu.txt > "$TMPDIR"/JOB_ID/out/cpu_tail.txt

rm "$TMPDIR"/JOB_ID/run/output
rm "$TMPDIR"/JOB_ID/out/timings.txt
rm "$TMPDIR"/JOB_ID/out/info.txt
rm "$TMPDIR"/JOB_ID/out/cpu.txt

cp -rf "$TMPDIR"/JOB_ID $HOME/EMRI/code/sims