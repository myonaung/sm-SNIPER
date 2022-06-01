data=/data/cephfs/punim1129/probe-cap_seq_data/sm-SNIPER/fastq

for file in ${data}/*.fastq;do
        outname=$(echo $file| rev | cut -d/ -f1 | rev )
        sample=$(echo $outname| cut -d. -f1)
        mkdir -p results
        mkdir -p results/PEPPER
        mkdir -p results/bam
        mkdir -p results/report
        mkdir -p results/${sample}
        mkdir -p results/${sample}/out
        mkdir -p results/${sample}/out/QC
        mkdir -p results/${sample}/out/confident_snv
        mkdir -p results/${sample}/out/tabular
        mkdir -p results/PEPPER/${sample}
        mkdir -p results/${sample}/out/text_file
        mkdir -p results/${sample}/out/snps_to_evaluate
        mkdir -p results/${sample}/out/snps_to_evaluate/longshot_pass
        mkdir -p results/${sample}/out/snps_to_evaluate/pepper_pass
        mkdir -p results/${sample}/out/snps_to_evaluate/svm_train_data
        mkdir -p results/${sample}/out/final_vcf_SNIPER
done