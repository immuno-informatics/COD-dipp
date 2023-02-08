python3.6 src/create_denovo_entries.py -input_denovo ../../tables/denovo_90ALC.tsv -input_blat ../../blat/data_gene3FT/denovo_humandb_filter.pslx -output output/denovo_90ALC_exons.tsv
python3.6 src/exons_denovo2mztab_lowmem.py output/denovo_90ALC_exons.tsv output/denovo_90ALC_exons.mzTAB

cd output
conda activate ~/plgg_iccvs/tools/ipip_conda_envs/687f7637

probam="$HOME/git_repositories/ipip/scripts/proBAMconvert/python/proBAM.py"
file="denovo_90ALC_exons.mzTAB"

python $probam \
    --name denovo_exons \
    --version 94 \
    --species homo_sapiens \
    --file $file \
    --conversion_mode proBAM_psm \
    --mismatches 1 \
    --sorting_order coordinate &

declare -i i
i=0
for file in *.mzTAB
do
    i+=1
    n=$(($i%6))
    python $probam \
        --name denovo_exons_$(basename $file .mzTAB) \
        --version 94 \
        --species homo_sapiens \
        --file $file \
        --conversion_mode proBAM_psm \
        --mismatches 1 \
        --sorting_order coordinate > $(basename $file .mzTAB).log 2>&1 &
    if [ $n -eq 0 ]
    then
        wait
    fi
done

conda activate ~/plgg_iccvs/tools/ipip_conda_envs/687f7637

probam="$HOME/git_repositories/ipip/scripts/proBAMconvert/python/proBAM.py"
file="open_search.mzTAB"

python $probam \
    --name open_search \
    --version 94 \
    --species homo_sapiens \
    --file $file \
    --conversion_mode proBAM_psm \
    --sorting_order coordinate > $(basename $file .mzTAB).log 2>&1 &

python ../../split_bam_by_sample.py -input_bam meta.sort.bam -input_psm ../../../tables/closed_search_PSMs.tsv -outdir bam_splits > bam_splits/log 2>&1 &
python ../../split_bam_by_sample.py -input_bam denono_exons_meta.bam  -input_psm ../../../tables/denovo_90ALC.tsv -outdir bam_splits > bam_splits/log 2>&1 &