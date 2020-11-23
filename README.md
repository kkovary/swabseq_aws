# swabseq_aws

## Run
docker build -t kylekovary:swabseq_aws .
docker run -d --rm --privileged -e PASSWORD=octant


## To do
- [ ] Get docker image on Amazon Elastic Container Service
- [ ] Get bs to install from Docker
- [ ] Connect output to cloud (GitHub, Google Drive, etc.)
- [ ] Output S2/S2_spike ratios in results file
- [x] Seach for all UDIs (no extra info about patients)
- [x] Download run from BaseSpace
