# Cluster ASD sequences at 30% identity (requires MMseqs2 on PATH).
# DATA_DIR should match DATA_DIR in config.py (default: <repo>/data).
DATA_DIR=${ALLOSTERIC_DATA_DIR:-data}

mmseqs easy-cluster $DATA_DIR/ASD_seqs.fasta $DATA_DIR/ASD_30_clusterRes tmp --min-seq-id 0.3
