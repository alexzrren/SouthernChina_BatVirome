rule urmap:
    input:
        fq1 = lambda wildcards: SAMPLE_INDEX[wildcards.sample_index][0],
        fq2 = lambda wildcards: SAMPLE_INDEX[wildcards.sample_index][1],
    output:
        fq1 = temp(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '01_urmap_r1.fq.gz')),
        fq2 = temp(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '01_urmap_r2.fq.gz')),
    params:
        fq1 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '.01_urmap_r1.fq.gz'),
        fq2 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '.01_urmap_r2.fq.gz'),
        sw = config['software']['urmap'],
        db = config['database']['urmap'],
        st = config['software']['samtools'],
        t  = config['resource']['urmap']['cpu'],
        sk = config['software']['seqkit'],
    priority: 2
    benchmark:
        os.path.join(config['output_dir'], 'logs', 'benchmark', '{sample_index}.01_urmap.bmk')
    shell:
        '''
        export OMP_NUM_THREADS={params.t}
        {params.sw} -threads {params.t} -map2 {input.fq1} -reverse {input.fq2} -ufi {params.db} -samout /dev/stdout |\
        {params.st} fastq -f 13 -1 {params.fq1} -2 {params.fq2} --threads {params.t} - 
        mv {params.fq1} {output.fq1}
        mv {params.fq2} {output.fq2}
        {params.sk} stat -j {params.t} {output.fq1} > {output.fq1}.stat
        {params.sk} stat -j {params.t} {output.fq2} > {output.fq2}.stat
        '''

rule fastp:
    input:
        fq1 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '01_urmap_r1.fq.gz'),
        fq2 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '01_urmap_r2.fq.gz'),
    output:
        fq1 = temp(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '02_fastp_r1.fq.gz')),
        fq2 = temp(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '02_fastp_r2.fq.gz')),
        html= os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '02_fastp.html'),
        json= os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '02_fastp.json'),
    params:
        fq1 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '.02_fastp_r1.fq.gz'),
        fq2 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '.02_fastp_r2.fq.gz'),
        sw = config['software']['fastp'],
        pm = config['params']['fastp'],
        t  = config['resource']['fastp']['cpu'],
        sk = config['software']['seqkit'],
    priority: 3
    benchmark:
        os.path.join(config['output_dir'], 'logs', 'benchmark', '{sample_index}.02_urmap.bmk')
    shell:
        '''
        {params.sw} {params.pm} --thread {params.t} -i {input.fq1} -o {params.fq1} -I {input.fq2} -O {params.fq2} --json {output.json} --html {output.html}
        mv {params.fq1} {output.fq1}
        mv {params.fq2} {output.fq2}
        {params.sk} stat -j {params.t} {output.fq1} > {output.fq1}.stat
        {params.sk} stat -j {params.t} {output.fq2} > {output.fq2}.stat
        '''

rule prinseq:
    input:
        fq1 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '02_fastp_r1.fq.gz'),
        fq2 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '02_fastp_r2.fq.gz'),
    output:
        fq1 = temp(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '03_prinseq_r1.fq.gz')),
        fq2 = temp(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '03_prinseq_r2.fq.gz')),
    params:
        fq1 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '.03_prinseq_r1.fq.gz'),
        fq2 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '.03_prinseq_r2.fq.gz'),
        sw = config['software']['prinseq'],
        pm = config['params']['prinseq'],
        t  = config['resource']['prinseq']['cpu'],
        sk = config['software']['seqkit'],
    priority: 4
    benchmark:
        os.path.join(config['output_dir'], 'logs', 'benchmark', '{sample_index}.03_prinseq.bmk')
    shell:
        '''
        {params.sw} {params.pm} -threads {params.t} -fastq {input.fq1} -fastq2 {input.fq2} -out_good {params.fq1} -out_good2 {params.fq2}
        mv {params.fq1} {output.fq1}
        mv {params.fq2} {output.fq2}
        {params.sk} stat -j {params.t} {output.fq1} > {output.fq1}.stat
        {params.sk} stat -j {params.t} {output.fq2} > {output.fq2}.stat
        '''

rule sortmeRNA:
    input:
        fq1 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '03_prinseq_r1.fq.gz'),
        fq2 = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '03_prinseq_r2.fq.gz'),
    output:
        fq1 = temp(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', 'cleanreads_r1.fq.gz')),
        fq2 = temp(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', 'cleanreads_r2.fq.gz')),
    params:
        sw = config['software']['sortmeRNA'],
        pm = config['params']['sortmeRNA'],
        t  = config['resource']['sortmeRNA']['cpu'],
        dbs= ' --ref '.join([config['database']['silva_arc_16s'], config['database']['silva_arc_23s'], config['database']['silva_bac_16s'], config['database']['silva_bac_23s'], config['database']['silva_euk_18s'], config['database']['silva_euk_28s']]),
        sk = config['software']['seqkit'],
        wd = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}'),
        pf1= os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '{sample_index}_norRNA'),
        pf2= os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', '{sample_index}_rRNA'),
    priority: 5
    benchmark:
        os.path.join(config['output_dir'], 'logs', 'benchmark', '{sample_index}.04_sortMeRNA.bmk')
    shell:
        '''
        if [ -e "{params.wd}/kvdb" ]; then
            rm -rf "{params.wd}/kvdb"
        fi
        if [ -e "{params.wd}/readb" ]; then
            rm -rf "{params.wd}/readb"
        fi
        {params.sw} {params.pm} --threads {params.t} --reads {input.fq1} --reads {input.fq2} --ref {params.dbs} \
            --workdir {params.wd} --aligned {params.pf2} --other {params.pf1}
        mv {params.pf1}_fwd.fq.gz {output.fq1}
        mv {params.pf1}_rev.fq.gz {output.fq2}
        {params.sk} stat -j {params.t} {output.fq1} > {output.fq1}.stat
        {params.sk} stat -j {params.t} {output.fq2} > {output.fq2}.stat
        if [ -e "{params.wd}/kvdb" ]; then
            rm -rf "{params.wd}/kvdb"
        fi
        if [ -e "{params.wd}/readb" ]; then
            rm -rf "{params.wd}/readb"
        fi
        if [ -e "{params.pf1}_fwd.fq.gz" ]; then
            rm -rf "{params.pf1}_fwd.fq.gz"
        fi
        if [ -e "{params.pf1}_rev.fq.gz" ]; then
            rm -rf "{params.pf1}_rev.fq.gz"
        fi
        '''

rule merge_cleanreads:
    input:
        expand(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', 'cleanreads_r1.fq.gz'), sample_index=SAMPLE_INDEX.keys()),
        expand(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'filter', '{sample_index}', 'cleanreads_r2.fq.gz'), sample_index=SAMPLE_INDEX.keys()),
    output:
        expand(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'merge', 'merge_reads', '{sample}', 'merged_cleanreads_r1.fq.gz'), sample=SAMPLE.keys()),
        expand(os.path.join(config['output_dir'], 'intermediate', 'upstream', 'merge', 'merge_reads', '{sample}', 'merged_cleanreads_r2.fq.gz'), sample=SAMPLE.keys()),
    params:
        sw = config['software']['merge_cleanreads'],
        wd = config['output_dir'],
    priority: 5
    benchmark:
        os.path.join(config['output_dir'], 'logs', 'benchmark', '05_merge_cleanreads.bmk')
    shell:
        '''
        {params.sw} --wd {params.wd}
        '''

for sample in SAMPLE.keys():
    SAMPLE[sample] = [os.path.join(config['output_dir'], 'intermediate', 'upstream', 'merge', 'merge_reads', f'{sample}', 'merged_cleanreads_r1.fq.gz'),
                      os.path.join(config['output_dir'], 'intermediate', 'upstream', 'merge', 'merge_reads', f'{sample}', 'merged_cleanreads_r2.fq.gz')]
