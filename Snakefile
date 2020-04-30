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
        config['variant_counts_file'],
#        config['func_scores_by_barcode_file']


rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        process_ccs=nb_markdown('process_ccs.ipynb'),
        build_variants=nb_markdown('build_variants.ipynb'),
        count_variants=nb_markdown('count_variants.ipynb'),
#        analyze_counts=nb_markdown('analyze_counts.ipynb')
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

rule analyze_counts:
    """Analyze variant counts and compute functional scores."""
    input:
        config['variant_counts_file']
    output:
        config['func_scores_by_barcode_file'],
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
