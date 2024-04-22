rule megahit:
    input:
        fq1 = lambda wildcards: SAMPLE[wildcards.sample][0],
        fq2 = lambda wildcards: SAMPLE[wildcards.sample][1],
    output:
        os.path.join(config['output_dir'], 'intermediate', 'upstream', 'assembly', '{sample}', '{sample}.contigs.fa')
    params:
        sw = config['software']['megahit'],
        pm = config['params']['megahit'],
        t  = config['resource']['megahit']['cpu'],
        sk = config['software']['seqkit'],
        wd = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'assembly', '{sample}'),
        wd_root = os.path.join(config['output_dir'], 'intermediate', 'upstream', 'assembly', '{sample}'),
        tmp= os.path.join(config['output_dir'], 'intermediate', 'upstream', 'assembly', '{sample}', 'intermediate_contigs'),
    priority: 6
    benchmark:
        os.path.join(config['output_dir'], 'logs', 'benchmark', '{sample}.06_megahit.bmk')
    shell:
        '''
        if [ -e "{params.wd}" ]; then
            rm -rf "{params.wd}"
        fi
        {params.sw} {params.pm} -t {params.t} --out-dir {params.wd} --out-prefix {wildcards.sample} \
            -1 {input.fq1} -2 {input.fq2}
        {params.sk} stat -j {params.t} {output} > {output}.stat
        if [ -e "{params.tmp}" ]; then
            rm -rf "{params.tmp}"
        fi
        '''