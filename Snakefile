"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import os.path
import textwrap

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: all,
            make_dag,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------
pacbio_runs = (pd.read_csv(config['pacbio_runs'], dtype = str)
               .assign(pacbioRun=lambda x: x['library'] + '_' + x['run'])
               )
assert len(pacbio_runs['pacbioRun'].unique()) == len(pacbio_runs['pacbioRun'])

# Rules -----------------------------------------------------------------------
rule all:
    """Final output of workflow."""
    input:
        os.path.join(config['summary_dir'], 'summary.md'),
        config['single_mut_effects_file']

rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        process_ccs=nb_markdown('process_ccs.ipynb'),
        build_variants=nb_markdown('build_variants.ipynb'),
        count_variants=nb_markdown('count_variants.ipynb'),
        analyze_counts=nb_markdown('analyze_counts.ipynb'),
        compute_Kd='results/summary/compute_binding_Kd.md',
        compute_meanF='results/summary/compute_expression_meanF.md',
        global_epistasis_binding=nb_markdown('global_epistasis_binding.ipynb'),
        global_epistasis_expression=nb_markdown('global_epistasis_expression.ipynb'),
        single_mut_effects='results/summary/single_mut_effects.md'
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:

            1. [Process PacBio CCSs]({path(input.process_ccs)}).

            2. [Build variants from CCSs]({path(input.build_variants)}).

            3. [Count variants by barcode]({path(input.count_variants)}).

            4. [QC analysis of sequencing counts]({path(input.analyze_counts)}).
            
            5. [Computation of ACE2-binding *K*<sub>D</sub>]({path(input.compute_Kd)}).
            
            6. [Computation of expression mean fluorescence]({path(input.compute_meanF)}).
            
            7. [Global epistasis decomposition of binding effects]({path(input.global_epistasis_binding)}).
            
            8. [Global epistasis decomposition of expression effects]({path(input.global_epistasis_expression)}).
            
            9. [Calculation of final single mutant effects on binding and expression]({path(input.single_mut_effects)}).

            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"

rule single_mut_effects:
    input:
        config['global_epistasis_binding_file'],
        config['global_epistasis_expr_file'],
        config['Titeseq_Kds_homologs_file'],
        config['expression_sortseq_homologs_file'],
    output:
        config['single_mut_effects_file'],
        config['homolog_effects_file'],
        md='results/summary/single_mut_effects.md',
        md_files = directory('results/summary/single_mut_effects_files')
    envmodules:
        'R/3.6.1-foss-2016b'
    params:
        nb='single_mut_effects.Rmd',
        md='single_mut_effects.md',
        md_files='single_mut_effects_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule global_epistasis_binding:
    input:
        config['Titeseq_Kds_file']
    output:
        config['global_epistasis_binding_file'],
        nb_markdown=nb_markdown('global_epistasis_binding.ipynb')
    params:
        nb='global_epistasis_binding.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule global_epistasis_expression:
    input:
        config['expression_sortseq_file']
    output:
        config['global_epistasis_expr_file'],
        nb_markdown=nb_markdown('global_epistasis_expression.ipynb')
    params:
        nb='global_epistasis_expression.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule compute_Titeseq_Kds:
    input:
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file'],
        config['Titeseq_Kds_homologs_file'],
        md='results/summary/compute_binding_Kd.md',
        md_files=directory('results/summary/compute_binding_Kd_files')
    envmodules:
        'R/3.6.1-foss-2016b'
    params:
        nb='compute_binding_Kd.Rmd',
        md='compute_binding_Kd.md',
        md_files='compute_binding_Kd_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule compute_expression_meanFs:
    input:
        config['variant_counts_file']
    output:
        config['expression_sortseq_file'],
        config['expression_sortseq_homologs_file'],
        md='results/summary/compute_expression_meanF.md',
        md_files=directory('results/summary/compute_expression_meanF_files')
    envmodules:
        'R/3.6.1-foss-2016b'
    params:
        nb='compute_expression_meanF.Rmd',
        md='compute_expression_meanF.md',
        md_files='compute_expression_meanF_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """


rule analyze_counts:
    """Analyze variant counts and compute functional scores."""
    input:
        config['variant_counts_file']
    output:
        nb_markdown=nb_markdown('analyze_counts.ipynb')
    params:
        nb='analyze_counts.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['codon_variant_table_file'],
        config['barcode_runs']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule build_variants:
    """Build variant table from processed CCSs."""
    input:
        config['processed_ccs_file']
    output:
        config['codon_variant_table_file'],
        nb_markdown=nb_markdown('build_variants.ipynb')
    params:
        nb='build_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule process_ccs:
    """Process the PacBio CCSs."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
        config['amplicons'],
        ([] if config['seqdata_source'] != 'HutchServer' else
         expand(os.path.join(config['ccs_dir'], "{pacbioRun}_report.txt"),
                pacbioRun=pacbio_runs['pacbioRun'])
         )
    output:
        config['processed_ccs_file'],
        nb_markdown=nb_markdown('process_ccs.ipynb')
    params:
        nb='process_ccs.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

if config['seqdata_source'] == 'HutchServer':

    rule build_ccs:
        """Run PacBio ``ccs`` program to build CCSs from subreads."""
        input:
            subreads=lambda wildcards: (pacbio_runs
                                        .set_index('pacbioRun')
                                        .at[wildcards.pacbioRun, 'subreads']
                                        )
        output:
            ccs_report=os.path.join(config['ccs_dir'], "{pacbioRun}_report.txt"),
            ccs_fastq=os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz")
        params:
            min_ccs_length=config['min_ccs_length'],
            max_ccs_length=config['max_ccs_length'],
            min_ccs_passes=config['min_ccs_passes'],
            min_ccs_accuracy=config['min_ccs_accuracy']
        threads: config['max_cpus']
        shell:
            """
            ccs \
                --min-length {params.min_ccs_length} \
                --max-length {params.max_ccs_length} \
                --min-passes {params.min_ccs_passes} \
                --min-rq {params.min_ccs_accuracy} \
                --report-file {output.ccs_report} \
                --num-threads {threads} \
                {input.subreads} \
                {output.ccs_fastq}
            """

elif config['seqdata_source'] == 'SRA':
    raise RuntimeError('getting sequence data from SRA not yet implemented')

else:
    raise ValueError(f"invalid `seqdata_source` {config['seqdata_source']}")
