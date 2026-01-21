SAMPLES = ["SRR2584403_1", "SRR2584404_1", "SRR2584405_1", "SRR2584857_1"]
GENOME = ["ecoli-rel606"]

rule make_vcf:
    input:
        expand("outputs/{sample}.x.{genome}.vcf",
               sample=SAMPLES, genome=GENOME),
        expand("outputs/{sample}.x.{genome}.vep.txt",
              sample=SAMPLES, genome=GENOME),
        "isec/0000.vcf"
  
rule uncompress_genome:
    input: "{genome}.fa.gz"
    output: "outputs/{genome}.fa"
    shell: """
        gunzip -c {input} > {output}
    """

rule map_reads:
    input:
        reads="{reads}.fastq.gz",
        ref="outputs/{genome}.fa"
    output: "outputs/{reads}.x.{genome}.sam"
    conda: "mapping"
    shell: """
        minimap2 -ax sr {input.ref} {input.reads} > {output}
    """

rule sam_to_bam:
    input: "outputs/{reads}.x.{genome}.sam",
    output: "outputs/{reads}.x.{genome}.bam",
    conda: "mapping"
    shell: """
        samtools view -b {input} > {output}
     """

rule sort_bam:
    input: "outputs/{reads}.x.{genome}.bam"
    output: "outputs/{reads}.x.{genome}.bam.sorted"
    conda: "mapping"
    shell: """
        samtools sort {input} > {output}
    """

rule index_bam:
    input: "outputs/{reads}.x.{genome}.bam.sorted"
    output: "outputs/{reads}.x.{genome}.bam.sorted.bai"
    conda: "mapping"
    shell: """
        samtools index {input}
    """

rule call_variants:
    input:
        ref="outputs/{genome}.fa",
        bam="outputs/{reads}.x.{genome}.bam.sorted",
        bai="outputs/{reads}.x.{genome}.bam.sorted.bai",
    output:
        pileup="outputs/{reads}.x.{genome}.pileup",
        bcf="outputs/{reads}.x.{genome}.bcf",
        vcf="outputs/{reads}.x.{genome}.vcf",
    conda: "mapping"
    shell: """
        bcftools mpileup -Ou -f {input.ref} {input.bam} > {output.pileup}
        bcftools call -mv -Ob {output.pileup} -o {output.bcf}
        bcftools view {output.bcf} > {output.vcf}
    """

rule tabix:
    input:
        gff="{filename}.gff.gz",
    output:
        tabix_idx='{filename}.gff.gz.tbi',
    conda: "mapping"
    shell: """
        tabix {input}
    """

rule predict_effects:
    input:
        fasta="{genome}.fa.gz",
        gff="{genome}.sorted.gff.gz",
        vcf="outputs/{reads}.x.{genome}.vcf",
        tabix_idx='ecoli-rel606.sorted.gff.gz.tbi',
    output:
        txt="outputs/{reads}.x.{genome}.vep.txt",
        html="outputs/{reads}.x.{genome}.vep.txt_summary.html",
        warn="outputs/{reads}.x.{genome}.vep.txt_warnings.txt",
    conda: "vep"
    shell: """
       vep --fasta {input.fasta} --gff {input.gff} -i {input.vcf} -o {output.txt}
    """

rule compress_vcfs:
    input:
        vcf="outputs/{reads}.x.{genome}.vcf"
    output:
        gz="isec/{reads}.x.{genome}.vcf.gz",
    conda: "mapping"
    shell: 
        """
        mkdir -p isec
        bgzip -kc {input} > {output}
        """

rule index_vcfs:
    input:
        gz="isec/{reads}.x.{genome}.vcf.gz"
    output:
        csi="isec/{reads}.x.{genome}.vcf.gz.csi"
    conda: "mapping"
    shell: 
        """
        bcftools index {input}
        """

rule compare_vcfs:
    input:
        gz=expand("isec/{sample}.x.{genome}.vcf.gz",
               sample=SAMPLES, genome=GENOME),
        csi=expand("isec/{sample}.x.{genome}.vcf.gz.csi",
              sample=SAMPLES, genome=GENOME),
    output:
        "isec/0000.vcf"
    conda: "mapping"
    shell:
        """
        bcftools isec -n +2 {input.gz} -p isec
        """