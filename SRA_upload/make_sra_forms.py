"""Script to upload FASTQ files to SRA.

Juhye Lee & Jesse Bloom, April 2019.
"""

import collections
import glob
import os
import tarfile

import yaml

import pandas as pd


def main():
    
    # output files 
    submissionform = 'submissionform_2019_Perth_serum_MAP.tsv'
    tar_dir = 'files_to_upload'
    os.makedirs(tar_dir, exist_ok=True)
    tar_filename = os.path.join(tar_dir, '2019_Perth_serum_MAP.tar')
    make_tar = True  # set to True to make the tarfile;

    # information for SRA about projects / samples
    bioprojectid = 'PRJNA530277'
    group_to_biosample = {  # assign serum groups to sample
                          'mock': 'SAMN11310373',
                          'plasmid': 'SAMN11310373',
                          'ferret': 'SAMN11310371',
                          'antibody': 'SAMN11310465',
                          'VIDD_sera': 'SAMN11310372',
                          'Hensley_sera': 'SAMN11310372',
                          'serum_mAb_spike': 'SAMN11341221',
                          }

    # read in all samples
    samplefile = '../data/sample_list.csv'
    samples = pd.read_csv(samplefile)

    # read in information sera
    serafile = '../data/serum_info.yaml'
    with open(serafile) as f:
        sera = (pd.DataFrame(yaml.safe_load(f))
               .transpose()
               .add_prefix('serum_')
               .rename_axis('serum')
               .reset_index()
               )

    # add sera information to samples
    missing_sera = set(samples['serum']) - set(sera['serum'])
    if missing_sera:
        raise ValueError(f"These serum not in {serafile}: {missing_sera}")
    samples = samples.merge(sera,
                            on='serum',
                            how='left',
                            validate='many_to_one'
                            )

    # add BioSample to samples
    missing_biosample = set(samples['serum_group']) - set(group_to_biosample)
    if missing_biosample:
        raise ValueError(f"No BioSample for serum groups: {missing_biosample}")
    samples = (samples
               .assign(biosample=lambda x: (x['serum_group']
                                            .map(group_to_biosample))
                       )
               )

    # Now build submission data frame, taking sample-specific columns
    # from `samples` data frame and assigning shared values to other columns.
    # Does not yet include files!
    submission = pd.DataFrame(collections.OrderedDict(
                    [('bioproject_accession', bioprojectid),
                     ('biosample_accession', samples['biosample']),
                     ('library_ID', samples['sample']),
                     ('title', samples['sample'] + ', ' + samples['serum_name']
                               + ' (' + samples['serum_description'] + ')'),
                     ('library_strategy', 'AMPLICON'),
                     ('library_source', 'VIRAL RNA'),
                     ('library_selection', 'PCR'),
                     ('library_layout', 'paired'),
                     ('platform', 'ILLUMINA'),
                     ('instrument_model', 'Illumina HiSeq 2500'),
                     ('design_description', '250-nt paired end reads of '
                                            'influenza HA PCR amplicons'),
                     ('filetype', 'fastq'),
                     ]))

    # Get all files, creating short names for each full path.
    file_d = collections.defaultdict(list)
    for tup in samples.itertuples():
        library_ID = tup.sample
        r1files = glob.glob(tup.R1)
        if len(r1files) == 0:
            raise ValueError(f"no R1 files matching {tup.R1}")
        r2files = [r1.replace('_R1_', '_R2_') for r1 in r1files]
        if not all(os.path.isfile(r2) for r2 in r2files):
            raise ValueError(f"missing R2 files: {r2files}")
        r1short = [f"{library_ID}_R1_{i + 1}.fastq.gz"
                   for i in range(len(r1files))]
        r2short = [f"{library_ID}_R2_{i + 1}.fastq.gz"
                   for i in range(len(r2files))]
        nfiles = len(r1files) + len(r2files)
        file_d['library_ID'] += [library_ID] * nfiles
        file_d['filelabel'] += [f"filename{i + 1}" for i in range(nfiles)]
        file_d['filename'] += r1short + r2short
        file_d['fullpath'] += r1files + r2files
    files = pd.DataFrame(file_d)
    if len(files['filename']) != len(files['filename'].unique()):
        raise ValueError('non-unique file names')

    # get files in wide data frame that can be merged into submission
    files_wide = (files
                  .pivot(values='filename',
                         columns='filelabel',
                         index='library_ID'
                         )
                  .rename(columns={'filename1': 'filename'})
                  .reset_index()
                  )
    submission = submission.merge(files_wide,
                                  on='library_ID',
                                  validate='one_to_one'
                                  )

    # write submission form
    print(f"\nWriting submission form to {submissionform}\n")
    submission.to_csv(submissionform, sep='\t', index=False)

    # make tar file
    if not make_tar:
        print('Not generating TAR file as `make_tar` is `False`.')
    else:
        print(f"Generating {tar_filename} (this may take a while)...")
        try:
            with tarfile.open(tar_filename, mode='w') as f:
                for i, tup in enumerate(files.itertuples()):
                    print(f"{i + 1} of {len(files)}: {tup.filename}...")
                    f.add(tup.fullpath, arcname=tup.filename)
            print(f"Added all files to {tar_filename}")
        except:
            if os.path.isfile(tar_filename):
                os.remove(tar_filename)
            raise

if __name__ == '__main__':
    main()

