keepvslist = open(config['keepvs_list']).read().splitlines()

vsname_clean = { vsname: vsname.split('/')[0] for vsname in keepvslist}

bampath = dict([tuple(line.split()) for line in open(config['quantbam_table']).read().splitlines()])


vs2bampath_dict = { key:dict() for key in vsname_clean.values() }
keep_bampath = {}
for line in open(config['quant_table']).read().splitlines():
	vs, sn = line.split()
	try:
		vs2bampath_dict[vsname_clean[vs]][sn] = bampath[sn]
		keep_bampath[sn] = bampath[sn]
	except:
		continue

vsname_clean = { v:k for k,v in vsname_clean.items() }


rule all:
	input:
		expand(os.path.join(config['wdir'], '7.distmatx_long',"{vsname}.csv"), vsname=vs2bampath_dict.keys())

rule bam_indexing:
	input:
		bam = lambda wildcards: keep_bampath[wildcards.sn]
	output:
		bam_symlink = temp(os.path.join(config['wdir'], "0.bamtmp/{sn}.bam")),
		bai = temp(os.path.join(config['wdir'], "0.bamtmp/{sn}.bam.bai"))
	shell:
		'''
		ln -s {input.bam} {output.bam_symlink}
		samtools index {output.bam_symlink}
		'''

rule fetch_bam:
	input:
		bam = rules.bam_indexing.output.bam_symlink,
		bai = rules.bam_indexing.output.bai
	output:
		os.path.join(config['wdir'], "1.fetchbam/{vsname}/{sn}.bam")
	params:
		contigid = lambda wildcards: vsname_clean[wildcards.vsname]
	shell:
		'''
		samtools view -bS {input.bam} {params.contigid} -o {output}
		'''

rule variant_calling:
	input:
		bamfiles = lambda wildcards: expand(os.path.join(config['wdir'], "1.fetchbam", f"{wildcards.vsname}", "{snlist}.bam"), snlist=vs2bampath_dict[wildcards.vsname].keys()),
		refgenomes = config['refgenomes']
	output:
		rglist = os.path.join(config['wdir'], '2.rglist', '{vsname}.rglist'),
		refgenome = os.path.join(config['wdir'], '2.rglist', '{vsname}.fna'),
		vcf = os.path.join(config['wdir'], '3.variant_calling',"{vsname}.vcf")
	params:
		contigid = lambda wildcards: vsname_clean[wildcards.vsname]
	shell:
		'''
		ls {input.bamfiles} > {output.rglist}
		seqkit grep -p {params.contigid} {input.refgenomes} > {output.refgenome}
		freebayes -f {output.refgenome} -L {output.rglist} -p 1 -0 > {output.vcf}
		'''

rule variant_qc:
	input:
		rules.variant_calling.output.vcf
	output:
		vcf_afterindvqc = os.path.join(config['wdir'], '4.vcf_indvQC', "{vsname}.vcf"),
		vcf_aftersiteqc = os.path.join(config['wdir'], '5.vcf_siteQC', "{vsname}.vcf"),
	params:
		remove_indvlist = os.path.join(config['wdir'], '4.vcf_indvQC', ".{vsname}/lowDP.indv"),
		tmpdir = os.path.join(config['wdir'], '4.vcf_indvQC', '.{vsname}')
	shell:
		'''
		mkdir -p {params.tmpdir}
		cd {params.tmpdir}
		vcftools --vcf {input} --missing-indv
		awk '$5 > 0.5' out.imiss | cut -f1 > {params.remove_indvlist}
		cd ..
		vcftools --vcf {input} --remove {params.remove_indvlist} --remove-filtered-all --recode --recode-INFO-all --stdout | vcfsnps  > {output.vcf_afterindvqc}
		vcftools --vcf {output.vcf_afterindvqc} --max-missing 0.75 --max-alleles 2 --minQ 30 --minDP 1 --max-alleles 2 --remove-filtered-all --recode --recode-INFO-all --stdout | vcfsnps -c > {output.vcf_aftersiteqc}
		
		if [ -e {params.tmpdir} ]; then
			rm -rf {params.tmpdir}
		fi
		'''

rule vcf2msa:
	input:
		rules.variant_qc.output.vcf_aftersiteqc
	output:
		tab = temp(os.path.join(config['wdir'], '6.consensus_SNPseq', '.{vsname}_popgenome.tab')),
		msa = os.path.join(config['wdir'], '6.consensus_SNPseq',"{vsname}.afa")
	params:
		vcf_tab2fas = '/jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Pipeline/viralstrain_phylo/scripts/vcf_tab_to_fasta_alignment.pl'
	shell:
		'''
		cat {input} | vcf-to-tab > {output.tab}
		{params.vcf_tab2fas} -i {output.tab} > {output.msa}
		'''

rule msa2distmatx:
	input:
		rules.vcf2msa.output.msa
	output:
		os.path.join(config['wdir'], '7.distmatx_long', "{vsname}.csv")
	params:
		msa2distmatx = '/jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Pipeline/viralstrain_phylo/scripts/msa2distmatx.py'
	shell:
		'''
		python {params.msa2distmatx} {input} {output}
		'''
