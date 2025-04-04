import pandas as pd
import os


# the problem arises when you try to generate the files from a master file in a bash loop. 
samples=pd.read_csv(config["sample_file"], header=None, sep=" ", quotechar="'",dtype=str)
samples.set_index(samples.iloc[:,0], inplace=True)

sampleIDs=samples.iloc[:,0].to_list()

reads=samples.iloc[0,1].split(",")

tissue=samples.iloc[0,2]
if '"' in tissue:
    tissue=tissue.strip('"')
print(tissue)

output_directory=config["root_dir"]
dir_to_samples = [output_directory]

qc_sum_output = list()
qc_sum_output.extend(expand(os.path.join("{dir}","{sampleID}.qcsum.txt"),zip, dir=dir_to_samples,sampleID=sampleIDs))

#check the present of hpo_term file
if any(samples.iloc[:,3].notna()):
    hpofiles=samples.iloc[:,3]
    hpo_file = hpofiles[hpofiles.notna()].unique()[0]
else:
    hpo_file = "none"

rule all:
    input:
        qc_sum=qc_sum_output
    params:
        outdir=output_directory,
        sampleID=sampleIDs[0],
        multiqc_prefix=sampleIDs[0]+'multiqc'
    output:
        multiqc_html="/".join([output_directory, sampleIDs[0]+"multiqc.html"])
    shell:
        """
        cd {params.outdir}

        mkdir err_out_folder
        find . -name '*.err' -exec mv -t ./err_out_folder/ {{}} +
        find . -name '*.out' -exec mv -t ./err_out_folder/ {{}} +

        multiqc {params.outdir} -o {params.outdir} -n {params.multiqc_prefix}
        """


rule trimming:
    input:
        r1=reads[0],
        r2=reads[1]
    output:
        paired1=temp(os.path.join('{dir}',"{sampleID}_R1_trim_paired.fq.gz")),
        paired2=temp(os.path.join('{dir}',"{sampleID}_R2_trim_paired.fq.gz")),
        unpaired1=os.path.join('{dir}',"{sampleID}_unpairR1.fq.gz"),
        unpaired2=os.path.join('{dir}',"{sampleID}_unpairR2.fq.gz"),
        failed_read=os.path.join('{dir}',"{sampleID}_failed.fq.gz"),
        json=os.path.join('{dir}',"{sampleID}fastp.json"),
        html=os.path.join('{dir}',"{sampleID}fastp.html")
    params:
        outdir='{dir}',
        sampleID='{sampleID}',
        fastp=config["executables"]["fastp"]
    shell:
        """
        {params.fastp} -w 10 -i {input.r1} -I {input.r2} -o {output.paired1} -O {output.paired2} \
        --unpaired1 {output.unpaired1} --unpaired2 {output.unpaired2} --failed_out {output.failed_read} \
        -h {output.html} -j {output.json}
        """

rule star:
    input:
        star_input1 = rules.trimming.output.paired1,
        star_input2 = rules.trimming.output.paired2
    output:
        tx_align=temp(os.path.join('{dir}',"{sampleID}Aligned.toTranscriptome.out.bam")),
        genome_align=temp(os.path.join('{dir}',"{sampleID}Aligned.out.bam")),
        junc_file=os.path.join('{dir}',"{sampleID}SJ.out.tab"),
        log_out=temp(os.path.join('{dir}',"{sampleID}Log.out")),
        log_final=os.path.join('{dir}',"{sampleID}Log.final.out"),
        progress_out=temp(os.path.join('{dir}',"{sampleID}Log.progress.out")),
        star_genome=temp(directory(os.path.join('{dir}',"{sampleID}_STARgenome"))),
        star_tmp=temp(directory(os.path.join('{dir}',"{sampleID}_STARtmp"))),
        star_pass=temp(directory(os.path.join('{dir}',"{sampleID}_STARpass1")))
    params:    
        index=config["resources"]["star_index"],
        gtf=config["resources"]["gtf"],
        STAR_modulefile=config["executables"]["STAR_modulefile"],
        prefix=os.path.join('{dir}',"{sampleID}"),
        outdir='{dir}',
        sampleID='{sampleID}'
    threads: 15
    shell:
        """
        module load {params.STAR_modulefile}; STAR --runMode alignReads \
            --runThreadN 15 \
            --readFilesCommand zcat \
            --readFilesIn {input.star_input1} {input.star_input2}\
            --genomeDir {params.index} \
            --outFileNamePrefix {params.prefix} \
            --twopassMode Basic \
            --sjdbGTFfile {params.gtf} \
            --outFilterType BySJout \
            --limitSjdbInsertNsj 1200000 \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs None \
            --alignSoftClipAtReferenceEnds Yes \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype BAM Unsorted \
            --outSAMattrRGline ID:{wildcards.sampleID}  SM:{wildcards.sampleID}  PL:ILLUMINA \
            --outSAMattributes All \
            --outSAMunmapped Within \
            --outSAMprimaryFlag AllBestScore \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --chimOutType Junctions \
            --chimMainSegmentMultNmax 1 \
            --genomeLoad NoSharedMemory
        """

rule sort:
    input:
        genome_align=rules.star.output.genome_align
    output:
        genome_align_sorted = temp(os.path.join('{dir}',"{sampleID}Aligned.sorted.bam"))
    params:
        outdir='{dir}',
        sampleID='{sampleID}'
    threads: 8
    shell:
        """
        module load samtools/1.9

        samtools sort -@ 8 -m 6G -o {output.genome_align_sorted} {input.genome_align}
        """

rule mark_duplicates:
    input:
        rules.sort.output.genome_align_sorted
    output:
        bam=os.path.join('{dir}',"{sampleID}.mdup.bam"),
        metrics=os.path.join('{dir}',"{sampleID}.duplicate.metrics")
    params:
        picard=config["executables"]["picard"],
        outdir='{dir}',
        sampleID='{sampleID}'
    threads: 10
    shell:
        """
        module load java/1.8.0_161
        module load samtools/1.9

        java -Xmx40g -jar {params.picard} MarkDuplicates \
        I={input} O={output.bam} M={output.metrics} \
        ASSUME_SORT_ORDER=coordinate 

        samtools index {output.bam}
        """

rule rsem:
    input:
        rules.star.output.tx_align
    output:
        gene=os.path.join('{dir}',"{sampleID}.genes.results"),
        isoform=os.path.join('{dir}',"{sampleID}.isoforms.results"),
        stat=temp(directory(os.path.join('{dir}',"{sampleID}.stat")))
    params:
        forward_prob="0",
        index=config["resources"]["rsem_index"],
        max_len="1000",
        outdir='{dir}',
        sampleID='{sampleID}',
        RSEM_modulefile=config["executables"]["RSEM_modulefile"]
    threads: 15
    shell:
        """
        module load {params.RSEM_modulefile}

        rsem-calculate-expression \
            --bam \
            --num-threads 15 \
            --fragment-length-max {params.max_len} \
            --no-bam-output \
            --paired-end \
            --forward-prob {params.forward_prob} \
            {input} \
            {params.index} {params.outdir}/{wildcards.sampleID}
        """

rule junctions:
    input:
        junctions=rules.star.output.junc_file,
        mdup_bam=rules.mark_duplicates.output.bam,
        database=ancient(config["resources"]["gtex_db"]),
        txdb=config["resources"]["tx_db"],
        annotations=config["resources"]["junctions"],
        refseq_db=config["resources"]["refseq_txdb"]
    params:
        tissue=tissue,
        hpo_file = hpo_file,
        outlier_frac=config["script_parameters"]["junction_outlier_frac"],
        outlier_mean_ratio=config["script_parameters"]["junction_outlier_mean_ratio"],
        novel_ratio=config["script_parameters"]["junction_novel_ratio"],
        novel_frac=config["script_parameters"]["junction_novel_frac"],
        z_cutoff=config["script_parameters"]["junction_zscore_cutoff"],
        missing_ratio=config["script_parameters"]["junction_missing_ratio"],
        missing_gtex_ratio=config["script_parameters"]["junction_missing_gtex_ratio"],
        missing_gtex_frac=config["script_parameters"]["junction_missing_gtex_frac"],
        junction_cutoff=config["script_parameters"]["junction_cutoff"],
        outdir='{dir}',
        sampleID='{sampleID}'
    output:
        filtered = os.path.join('{dir}',"{sampleID}.junctions.tsv"),
        alljuncs = os.path.join('{dir}',"{sampleID}.junctions.all.tsv.gz")
    threads: 10
    shell:
        """
        module load R/4.0.3
        module load bedtools/2.29.2

        echo "{params.tissue}"

        Rscript scripts/junction_outlier.R -d {input.database} -t "{params.tissue}" -o {output.filtered} \
            -g {input.txdb} -a {input.annotations} -j {input.junctions} -s {wildcards.sampleID} --HPO_terms {params.hpo_file} --refseq_txdb {input.refseq_db} \
            --outlier_frac {params.outlier_frac} --outlier_mean_ratio {params.outlier_mean_ratio} --novel_ratio {params.novel_ratio} \
            --novel_frac {params.novel_frac} --z_cutoff {params.z_cutoff} --missing_ratio {params.missing_ratio} \
            --missing_gtex_ratio {params.missing_gtex_ratio} --missing_gtex_frac {params.missing_gtex_frac} --junction_cutoff {params.junction_cutoff}
        """

rule rnaseq_metrics:
    input:
        rules.mark_duplicates.output.bam
    output:
        metrics=os.path.join('{dir}',"{sampleID}.rnaseq_metrics.txt")
    params:
        ribo_int=config["resources"]["ribo_int"],
        refflat=config["resources"]["refflat"],
        strand="SECOND_READ_TRANSCRIPTION_STRAND",
        AS="true",
        picard=config["executables"]["picard"],
        outdir='{dir}',
        sampleID='{sampleID}'
    threads: 10
    shell:
        """
        module load java/1.8.0_161

        java -Xmx4096m -jar {params.picard} CollectRnaSeqMetrics \
        REF_FLAT={params.refflat} \
        RIBOSOMAL_INTERVALS={params.ribo_int} \
        STRAND_SPECIFICITY= {params.strand} \
        I={input} \
        O={output.metrics} AS={params.AS}
        """

rule rnaseqc:
    input:
        bam = rules.mark_duplicates.output.bam,
        collapsed_gtf = config["resources"]["collapsed_gtf"]
    output:
        os.path.join('{dir}',"{sampleID}.metrics.tsv")
    params:
        rnaseqc = config["executables"]["rnaseqc"],
        outdir='{dir}',
        sampleID='{sampleID}'
    threads: 10
    shell:
        """
        {params.rnaseqc} {input.collapsed_gtf} {input.bam} {params.outdir} -s {wildcards.sampleID}
        rm -f {params.outdir}/*gct
        """

rule insert_size_metrics:
    input:
        rules.mark_duplicates.output.bam
    output:
        metrics=os.path.join('{dir}',"{sampleID}.insert_metrics.txt"),
        plot=temp(os.path.join('{dir}',"{sampleID}.insert_metrics.pdf")),
    params:
        picard=config["executables"]["picard"]
    threads: 10
    shell:
        """
        module load java/1.8.0_161
        module load R/3.5.1

        java -Xmx4096m -jar {params.picard} CollectInsertSizeMetrics \
        I={input} O={output.metrics} H={output.plot}
        """

rule ercc_sirv:
    input:
        gene=rules.rsem.output.gene,
        iso=rules.rsem.output.isoform,
        ercc_conc=config["resources"]["spikeins"],
        python_modulefile=config["executables"]["python_modulefile"]
    output:
        ercc=os.path.join('{dir}',"{sampleID}.ercc.txt"),
        sirv=os.path.join('{dir}',"{sampleID}.sirv.txt")
    threads: 10
    shell:
        """
        module load {input.python_modulefile}

        source /hpf/largeprojects/ccmbio/yliang/clinical_pipeline/python_venv/bin/activate

        python3 scripts/spikeins.py --gene {input.gene} --iso {input.iso} --actual {input.ercc_conc} \
            --samplename {wildcards.sampleID} --analysis_path {wildcards.dir}
        """

rule fastqc:
    input:
        #reads=lambda wildcards: samples.loc[wildcards.sampleID, 1].split(',')
        read_input1 = rules.trimming.output.paired1,
        read_input2 = rules.trimming.output.paired2
    output:
        report1=os.path.join('{dir}', "{sampleID}_1_fastqc/fastqc_data.txt"),
        report2=os.path.join('{dir}', "{sampleID}_2_fastqc/fastqc_data.txt")
    threads: 2
    shell:
        """
        module load fastqc/0.11.5

        fastqc -t 2 -q --extract -o {wildcards.dir} {input.read_input1} {input.read_input2}

        current_dirname1={wildcards.dir}/$(basename {input.read_input1} | cut -d '.' -f1)_fastqc
        target_dirname1=$(dirname {output.report1})

        current_dirname2={wildcards.dir}/$(basename {input.read_input2} | cut -d '.' -f1)_fastqc
        target_dirname2=$(dirname {output.report2})

        if [[ "$current_dirname1" != "$target_dirname1" ]]
        then
          rm -r $target_dirname1
          rm -r $target_dirname2
          mv $current_dirname1 $target_dirname1
          mv $current_dirname2 $target_dirname2
        fi
        """

rule cal_q30:
    input:
        read1=lambda wildcards: samples.loc[wildcards.sampleID, 1].split(',')[0],
        read2=lambda wildcards: samples.loc[wildcards.sampleID, 1].split(',')[1]
    output:
        qchist = os.path.join('{dir}',"{sampleID}.qchist.txt")
    params:
        outdir='{dir}',
        sampleID='{sampleID}'
    threads: 10
    shell:
        """
        module load bbmap/37.33; reformat.sh -Xmx1g in1={input.read1} in2={input.read2} qchist={output.qchist}
        """

rule get_cov:
    input:
        bam = rules.mark_duplicates.output.bam,
        genome_index = config["resources"]["fai"],
        bed = config["resources"]["exons_bed"]
    output:
        cov_file = os.path.join('{dir}',"{sampleID}.coverage.stat")
    params:
        outdir='{dir}',
        sampleID='{sampleID}'
    threads: 10
    shell:
        """
        module load samtools/1.9

        scripts/get_rnaseq_cov.sh {input.bam} {input.genome_index} {input.bed}
        """

rule alignment_metrics:
    input:
        rules.mark_duplicates.output.bam
    output: os.path.join('{dir}',"{sampleID}.alignment_metrics.txt")
    params:
        picard=config["executables"]["picard"],
        fasta=config["resources"]["fasta"],
        outdir='{dir}',
        sampleID='{sampleID}'
    threads: 10
    shell:
        """
        module load java/1.8.0_161

        java -Xmx4096m -jar {params.picard} CollectAlignmentSummaryMetrics \
        I={input} R={params.fasta} O={output}        
        """

rule samtools_idx:
    input:
        rules.mark_duplicates.output.bam
    output:
        os.path.join('{dir}',"{sampleID}.idxstats")
    threads: 10
    shell:
        """
        module load samtools/1.9

        samtools idxstats {input} > {output}
        """


rule sum_metrics:
    input:
        star_file = rules.star.output.log_final,
        junction_analysis_files = rules.junctions.output.filtered,
        bam_file = rules.mark_duplicates.output.bam,
        picard_file = rules.rnaseq_metrics.output.metrics,
        duplicate_file = rules.mark_duplicates.output.metrics,
        insert_file = rules.insert_size_metrics.output.metrics,
        alignment_file = rules.alignment_metrics.output,
        qchist_file = rules.cal_q30.output.qchist,
        coverage_file = rules.get_cov.output.cov_file,
        ercc_file = rules.ercc_sirv.output.ercc,
        sirv_file = rules.ercc_sirv.output.sirv,
        junc_file = rules.star.output.junc_file,
        rnaseq_file = rules.rnaseqc.output,
        expr_file = rules.rsem.output.gene,
        samtool_idxstat= rules.samtools_idx.output,
        fastqc1=rules.fastqc.output.report1,
        fastqc2=rules.fastqc.output.report2
    params:
        out_suffix=config["script_parameters"]["sum_qc_out_suffix"],
        seq_len=config["script_parameters"]["sum_qc_seq_len"],
        tx_len=config["script_parameters"]["sum_qc_tx_len"],
        outdir='{dir}',
        sampleID='{sampleID}'
    output:
        qc_sum=os.path.join('{dir}',"{sampleID}.qcsum.txt")
    threads: 10
    shell:
        """
        module load R/3.5.1
        module load samtools/1.9
        
         Rscript scripts/sum_qc_metrics.r -s {wildcards.sampleID} -d {params.outdir} -o {params.out_suffix} \
            --star_file {input.star_file} --bam_file {input.bam_file} --picard_file {input.picard_file} \
            --duplicate_file {input.duplicate_file} --insert_file {input.insert_file} --alignment_file {input.alignment_file} \
            --qchist_file {input.qchist_file} --coverage_file {input.coverage_file} --expr_file {input.expr_file} \
            --junc_file {input.junc_file} --ercc_file {input.ercc_file} --sirv_file {input.sirv_file} \
            --rnaseq_file {input.rnaseq_file} --seq_len {params.seq_len} --tx_len {params.tx_len}

        """
