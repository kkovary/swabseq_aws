# swabseq_aws

## Issues
* Docker is currently not working, issue with automated authentication of BaseSpace CLI in order to download data.
* Because of the issue described above, the docker file has not been updated to work with current code.


## Usage
 * `Rscript countAmpliconsAWS.R --basespaceID [ID for run] --threads [number of threads for running bcl2fastq]`
 * The basespaceID is used to identify the run on BaseSpace and then download the raw data which is then demultiplexed with bcl2fastq and then analyzed, where a PDF of run info and results is generated, along with a csv file with the unique DNA barcodes for each sample, the location of that sample on 96 and 384 well plates, the number of counts for the targeted amplicons, and the classification of the sample (COVID positive, COVID negative, or inconclusive/failed sample).

