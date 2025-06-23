# sh file that, for each bbnet in the directory, calls "train.sh in=ncbi_100_6f.tsv netin=***.bbnet evaluate"

for net in *.bbnet; do
    echo "Training and evaluating 'ncbi_100_6f' with network $net"
    ./bbmap/train.sh in=ncbi_100_6f.tsv netin="$net" evaluate
done