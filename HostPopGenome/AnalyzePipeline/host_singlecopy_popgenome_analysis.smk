with open(config['samplelist'], 'r') as fd:
	SAMPLE_INDEX = { x[0]:x[1:] for x in list(map(lambda x: x.split(), fd.read().splitlines()))}

with open(config['reference_fasta'], 'r') as fd:
    REFCDS = { str(i) : seqid for i, seqid in enumerate(list(map(lambda x: x.split()[0].replace('>', ''), list(filter(lambda x: x.startswith('>'), fd.read().splitlines()))))) }

#print(REFCDS)

from pathlib import Path

Path('logs').mkdir(exist_ok=True)

rule all:
    input:
        os.path.join(config['outdir'], '3.vcf_merge_QC', '3.4.merged_allcds_rmsingleton.vcf'),
        expand(os.path.join(config['outdir'], '4.PCA', 'singlecopy_cds.{suffixs}'), suffixs=['eigenval', 'eigenvec', 'ped', 'bed', 'fam'])


localrules: generate_rglist, plink_pca

rule bowtie2_mapping:
    input:
        fq1 = lambda wildcards: SAMPLE_INDEX[wildcards.sample_index][0],
        fq2 = lambda wildcards: SAMPLE_INDEX[wildcards.sample_index][1],
    params:
        reference_index = config['reference_index'],
        bowtie2 = config['software']['bowtie2'],
        rgid = lambda wildcards: wildcards.sample_index
    output:
        bam = os.path.join(config['outdir'], "1.sortedbam", "{sample_index}_sorted.bam"),
        bai = os.path.join(config['outdir'], "1.sortedbam", "{sample_index}_sorted.bam.bai")
    shell:
        """
        {params.bowtie2} -x {params.reference_index} -1 {input.fq1} -2 {input.fq2} --very-sensitive-local -p 4 | samtools addreplacerg -r "ID:{params.rgid}" -r "SM:{params.rgid}" - | samtools view -bS -F 4 | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule generate_rglist:
    input:
        expand(os.path.join(config['outdir'], '1.sortedbam', "{sample_index}_sorted.bam"), sample_index=SAMPLE_INDEX.keys())
    output:
        os.path.join(config['outdir'], 'allbamfiles.rglist')
    shell:
        "ls {input} > {output}"

rule freebayes_eachcds:
    input:
        reference_fasta = config['reference_fasta'],
        rgfile = rules.generate_rglist.output
    params:
        cds = lambda wildcards: REFCDS[wildcards.region_id]
    output:
        vcf = os.path.join(config['outdir'], '2.cds_variantcalling', 'cds_{region_id}.vcf'),
        milestone = os.path.join(config['outdir'], '.2.cds_variantcalling_milestone', 'cds_{region_id}.done')
    shell:
        '''
        freebayes -f {input.reference_fasta} -L {input.rgfile} -p 2 -r '{params.cds}' -0 > {output.vcf}
        touch {output.milestone}
        '''

rule merge_allvcf:
    input:
        expand(os.path.join(config['outdir'], '2.cds_variantcalling', 'cds_{region_id}.vcf'), region_id=REFCDS.keys())
    output:
        os.path.join(config['outdir'], '3.vcf_merge_QC', '3.1.merged_allcds.vcf')
    params:
        max_cdsid = len(REFCDS)-1,
        vcfpath = os.path.join(config['outdir'], '2.cds_variantcalling')
    shell:
        '''
        grep '#' {params.vcfpath}/cds_0.vcf > {output}
        for i in {{0..{params.max_cdsid}}}
        do
            grep -v '#' {params.vcfpath}/cds_$i.vcf >> {output}
        done
        '''

rule variant_qc:
    input:
        rules.merge_allvcf.output
    output:
        vcf_aftersiteqc = os.path.join(config['outdir'], '3.vcf_merge_QC', '3.2.merged_allcds_siteQC.vcf'),
        vcf_afterindvqc = os.path.join(config['outdir'], '3.vcf_merge_QC', '3.3.merged_allcds_indvQC.vcf'),
        vcf_rmsingleton = os.path.join(config['outdir'], '3.vcf_merge_QC', '3.4.merged_allcds_rmsingleton.vcf')
    params:
        remove_indvlist = os.path.join(config['outdir'], '.vcfQC_tmp/lowDP.indv'),
        tmpdir = os.path.join(config['outdir'], '.vcfQC_tmp'),
    shell:
        '''
        vcftools --vcf {input} --max-missing 0.5 --max-alleles 2 --minQ 30 --minDP 1 --max-alleles 2 --remove-filtered-all --recode --recode-INFO-all --stdout | vcfsnps -c > {output.vcf_aftersiteqc}

        mkdir -p {params.tmpdir}
        cd {params.tmpdir}
        vcftools --vcf {output.vcf_aftersiteqc} --missing-indv
        awk '$5 > 0.5' out.imiss | cut -f1 > {params.remove_indvlist}
        cd ..
        vcftools --vcf {output.vcf_aftersiteqc} --remove {params.remove_indvlist} --remove-filtered-all --recode --recode-INFO-all --stdout | vcfsnps  > {output.vcf_afterindvqc}
        if [ -e {params.tmpdir} ]; then
            rm -rf {params.tmpdir}
        fi

        vcftools --vcf {output.vcf_afterindvqc} --mac 3 --recode --recode-INFO-all --stdout > {output.vcf_rmsingleton}
        '''

rule plink_pca:
    input:
        rules.variant_qc.output.vcf_rmsingleton
    output:
        expand(os.path.join(config['outdir'], '4.PCA', 'singlecopy_cds.{suffixs}'), suffixs=['eigenval', 'eigenvec', 'ped', 'bed', 'fam'])
    params:
        output_prefix = os.path.join(config['outdir'], '4.PCA', 'singlecopy_cds')
    shell:
        '''
        plink --allow-extra-chr --vcf {input} --recode --out {params.output_prefix}
        plink --allow-extra-chr --file {params.output_prefix} --noweb --make-bed --out {params.output_prefix}
        plink --allow-extra-chr --threads 20 --bfile {params.output_prefix} --pca 20 --out {params.output_prefix}
        '''
