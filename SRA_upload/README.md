# Uploading FASTQ files to the SRA
The script instructions and for how the raw sequencing files (FASTQ files) were uploaded to the NIH [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra):

1. Go to the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) and manually create a *BioProject* and appropriate *BioSamples*.
   These are then entered into the Python script [make_sra_forms.py](make_sra_forms.py) near the top to specify the *BioProject* and *BioSamples* for the submission form that script automatically creates.

2. Run [make_sra_forms.py](make_sra_forms.py), a Python script that reads in the [../data/serum_info.yaml](../data/serum_info.yaml) and [../data/sample_list.csv](../data/sample_list.csv) files and uses them to generate the files needed for submission.
   Specifically, it generates the following two files in the subdirectory [./files_to_upload](files_to_upload):

     - the submission form [submissionform_2019_Perth_serum_MAP.tsv](submissionform_2019_Perth_serum_MAP.tsv) 

     - the FASTQ `tar` file `./files_to_upload/2019_Perth_serum_MAP.tar`.

   The `tar` file takes a long time to generate, so you might first want to run the script with the `make_tar` option set to `False` to make sure it is working.

3. Go to the [SRA submission portal](https://submit.ncbi.nlm.nih.gov/subs/sra/) and login.
   You will see three options to preload data: the Aspera browser plugin upload, the Aspera Command-Line upload, and FTP upload. 
   The instructions here are for the **Aspera Command-Line**, which appears to be the fastest method.

4. Click on the **Aspera Command-Line** option.
   The instructions here give the following details for this particular upload:

        ascp -i <path/to/key_file> -QT -l100m -k1 -d <path/to/folder/containing files> subasp@upload.ncbi.nlm.nih.gov:uploads/jbloom_fhcrc.org_8LXpJe7P

   Our script put the `tar` in the directory `files_to_upload` as described above.
   Note that we do **not** put the `*.tsv` submission form in the `files_to_upload` directory as that is uploaded later as described below.
   On the Hutch server, we want to use the version of `ascp` at `/app/aspera-connect/3.7.5/bin/ascp`.
   The key for this version of `ascp` is at `/app/aspera-connect/3.7.5/etc/asperaweb_id_dsa.openssh`.
   So we run the following command:

        /app/aspera-connect/3.7.5/bin/ascp -i /app/aspera-connect/3.7.5/etc/asperaweb_id_dsa.openssh -QT -l100m -k1 -d files_to_upload subasp@upload.ncbi.nlm.nih.gov:uploads/jbloom_fhcrc.org_8LXpJe7P

   On the Hutch cluster, it is easiest to submit this command via `sbatch` so it does not time out if your window closes:

        sbatch --wrap="/app/aspera-connect/3.7.5/bin/ascp -i /app/aspera-connect/3.7.5/etc/asperaweb_id_dsa.openssh -QT -l100m -k1 -d files_to_upload subasp@upload.ncbi.nlm.nih.gov:uploads/jbloom_fhcrc.org_8LXpJe7P"

   If you submit the upload via `sbatch`, then just use `squeue -u <YOUR USERNAME>` to see when upload is done, or look in the `slurm-*` file, which will print a message when upload is complete.
   The upload will take a while; for this particular upload (122 GB) it took about 2.5 hours.

5. After the upload is complete, return back to [SRA submission portal](https://submit.ncbi.nlm.nih.gov/subs/sra/) and click the **New submission** button:

  a. Fill out your information in the first page (I used "Fred Hutchinson Cancer Research Center" in the Submitting organization field and "Division of Basic Sciences" in the Department field)

  b. For 'Do you want to create new BioProject?', select 'No', and enter the existing BioProject ID. Select 'No' for 'Do you want to create new BioSamples for this submission?'

  c. Upload the `.tsv` submission on the SRA metadata page

  d. On the Files page, select "FTP or Aspera preload".
     Click on the `Select preload folder` button, and a pop-up will show up with the directory you uploaded the `.tar` file into.
     Select this directory.
     If you uploaded by the Aspera Command-line, the directory will be signified with '(Aspera)'.
     If you uploaded by FTP, the directory will be signified with '(FTP)'.
     Note that you will need to wait at least 10-15 minutes after the file uploads for the folder to show up here.
     Check Autofinish submission.
     Click 'Continue'.
     You will get an warning message that the files are missing as well as an automated e-mail warning you that files are missing, so click 'Continue' again to extract the `fastq` files from the `.tar` file.
     This can take anywhere between several minutes to an hour.
     If any `FASTQ` files are missing, you will receive an error message listing the missing files.

  e. Once all of the files are processed, you will be redirected to a page that says your samples have been submitted and are awaiting processing. The submission portal should also now show `Processing` under Status.
  
  f. You're now done. It will just take a little while for the files to become available.
