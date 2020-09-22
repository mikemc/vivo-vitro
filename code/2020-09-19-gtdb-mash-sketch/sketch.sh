# To sketch with mash on the cluster:
cd ~/data/gtdb
sed -e 's#^#genomes/#' file-list.txt > mash-file-list.txt
sbatch --exclusive --job-name=gtdb-mash --wrap="
  mash sketch -k 31 -s 2000 -l mash-file-list.txt -o gtdb-v95-k31s2000 -p 16
"
