"""Run Jupyter notebook and create Markdown output.

Run with:

    python run_nb.py <nb> <md_output>
    
"""


import argparse
import os
import subprocess


def main():

    parser = argparse.ArgumentParser(
                description='Run Jupyter notebook, create Markdown output')
    parser.add_argument('nb', help='Name of Jupyter notebook.')
    parser.add_argument('md_output', help='Name of created Markdown output.')
    args = parser.parse_args()

    nb = args.nb
    if not (os.path.isfile(nb) and os.path.splitext(nb)[-1] == '.ipynb'):
        raise IOError(f"not a valid Jupyter notebook: {nb}")

    md_output = args.md_output
    if os.path.splitext(md_output)[-1] not in {'.md', '.markdown'}:
        raise IOError(f"not a valid Markdown file: {md_output}")

    if os.path.isfile(md_output):
        os.remove(md_output)

    subprocess.check_call(['jupyter', 'nbconvert',
                           '--to', 'notebook',
                           '--execute',
                           '--inplace',
                           '--ExecutePreprocessor.timeout=-1',
                           nb,
                           ])

    subprocess.check_call(['jupyter', 'nbconvert',
                           '--output-dir', os.path.dirname(md_output),
                           '--output', os.path.basename(md_output),
                           '--to', 'markdown',
                           nb,
                           ])


if __name__ == '__main__':
    main()
