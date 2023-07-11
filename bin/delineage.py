import argparse
import pandas as pd
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import Normalizer
from scipy.cluster.hierarchy import linkage, to_tree
from Bio import SeqIO

def clustering_contigs(coverage_df, tetra_df):

    # Calculate the pairwise distance matrix
    tetra_distances = pairwise_distances(tetra_df, metric="euclidean")

    # Apply MDS and transform the data to 2D
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0)
    tetra_2d = mds.fit_transform(tetra_distances)
    tetra_2d_df = pd.DataFrame(tetra_2d, index=tetra_df.index, columns=["MDS1", "MDS2"])

    # Apply the normalizer to the coverage data
    normalizer = Normalizer(norm="l1")
    coverage_normalized = normalizer.transform(coverage_df)
    coverage_normalized_df = pd.DataFrame(coverage_normalized, index=coverage_df.index, columns=coverage_df.columns)

    # Concatenate the coverage and tetra dataframes
    combined_df = pd.concat([coverage_normalized_df, tetra_2d_df], axis=1)

    # Perform hierarchical clustering
    linked = linkage(combined_df, method='ward', metric='euclidean')

    # Convert the linkage matrix to a tree
    tree = to_tree(linked)
    contigid_list = list(combined_df.index)

    return tree, contigid_list

def get_newick(node, newick, labels):
    if node.is_leaf():
        return "%s%s" % (labels[node.id], newick)
    else:
        if len(newick) > 0:
            newick = "):%s%s" % (node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), newick, labels)
        newick = get_newick(node.get_right(), ",%s" % newick, labels)
        newick = "(%s" % newick
        return newick
    
#Get all leaves from a ClusterNode object
def get_leaves_from_cluster(cluster):
    if cluster.is_leaf():
        return [cluster.id]
    else:
        return get_leaves_from_cluster(cluster.get_left()) + get_leaves_from_cluster(cluster.get_right())

# Function to check requirements for each group of contigs
def check_requirements(group, contig_lengths, hallmark_counts, single_copy_genes, NCLDV_assessment):
    lengths = [contig_lengths[id] for id in group]
    total_length = sum(lengths)

    group_hallmarks = hallmark_counts[hallmark_counts['ID'].isin(group)]
    has_hallmark_gene = (group_hallmarks[single_copy_genes] > 0).any(axis=1).any()

    group_scores = NCLDV_assessment[NCLDV_assessment['ID'].isin(group)]
    over_2_percent = (group_scores['NCLDV_score'] > 2).mean()
    return total_length, has_hallmark_gene, over_2_percent

def create_summary_df(group1_reqs, group2_reqs, bin_id):
    data = {
        "ID": [bin_id, bin_id],
        "group": ["group1", "group2"],
        "total_length": [group1_reqs[0], group2_reqs[0]],
        "has_scgs": [int(group1_reqs[1]), int(group2_reqs[1])],
        "NCLDV_score_over_2_percent": [group1_reqs[2], group2_reqs[2]]
    }

    df = pd.DataFrame(data)

    # Add a new column for whether all criteria are met
    df["all_criteria_met"] = (df["total_length"] > 40000) & (df["has_scgs"] == 1) & (df["NCLDV_score_over_2_percent"] > 0.5)
    df["all_criteria_met"] = df["all_criteria_met"].astype(int)

    return df

def process_bin(summary_df, sequences, group1, group2):
    if summary_df["all_criteria_met"].all():  # If all criteria are met for both groups
        seq_group = split_fasta(sequences, group1, group2)  # Split the fasta file
        return seq_group
    else:
        annotated_seq = annotate_fasta(sequences)
        return annotated_seq

def split_fasta(sequences, group1, group2):
    # Split the sequences into two groups
    group1_sequences = [seq for seq in sequences if seq.id in group1]
    group2_sequences = [seq for seq in sequences if seq.id in group2]
    return [group1_sequences, group2_sequences]


def annotate_fasta(sequences):
    # Add the annotation to each sequence description
    for seq in sequences:
        seq.description += " multiple SCG NCLDV bin"
    return [sequences]

def main():
    """Process command line arguments and run the script"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bin_id", required=True, type=str, help="The bin id")
    parser.add_argument("-f", "--fasta", required=True, type=str, help="Path to the bin contigs fasta file")
    parser.add_argument("-t", "--tetranucleotide", required=True, type=str, help="Path to the 4-mer counts TSV file (produced by cgat)")
    parser.add_argument('-d', '--depth', required=True, type=str, help="Path to the depth file (totalAvgDepth in 3rd column) in TSV format.")
    parser.add_argument('-m', '--hallmark_gene_summary', required=True, type=str, help="Path to the hallmark gene HMMER results summary per contig file in TSV format.")
    parser.add_argument("-n", "--ncldv_assessment", required=True, type=str, help="Path to the NCLDV score assessemt TSV file")
    parser.add_argument('-s', '--scgs', required=True, type=str, help='Comma-separated list of NCLDV single copy genes.')
    parser.add_argument("-o", "--output", required=True, type=str, help="Output directory for bins")
    parser.add_argument("-O", "--tree_output", required=True, type=str, help="Path to the tree output")
    parser.add_argument("-S", "--summary_output", required=True, type=str, help="Path to the summary output TSV file.")

    args = parser.parse_args()

    # Convert the single-copy genes to a list
    single_copy_genes = args.scgs.split(',')

    coverage_df = pd.read_csv(args.depth, sep="\t", header=0)
    coverage_df = coverage_df.iloc[:, [0] + list(range(3, len(coverage_df.columns), 2))]
    coverage_df.set_index("contigName", inplace=True)
    # Identify all-zero columns and drop all-zero columns from the DataFrame
    coverage_df = coverage_df.loc[:, ~(coverage_df == 0).all()]

    # Load the 4-mer counts TSV file, transpose it, and set contig_id as index
    tetra_df = pd.read_csv(args.tetranucleotide, sep="\t", header=0).T
    new_header = tetra_df.iloc[0]
    tetra_df = tetra_df[1:]
    tetra_df.columns = new_header
    tetra_df.index.name = "contigName"
    # Identify all-zero columns and drop all-zero columns from the DataFrame
    tetra_df = tetra_df.loc[:, ~(tetra_df == 0).all()]

    # get clustering result in tree format
    tree, contigid_list = clustering_contigs(coverage_df, tetra_df)

    # Convert the tree to Newick format
    newick = get_newick(tree, "", contigid_list)

    # Write the Newick string to a file
    with open(args.tree_output, 'w') as f:
        f.write(newick)

    # Extract the contig IDs from the two major branches of the tree
    group1 = [contigid_list[node] for node in get_leaves_from_cluster(tree.get_left())]
    group2 = [contigid_list[node] for node in get_leaves_from_cluster(tree.get_right())]

    # Load the fasta file
    sequences = [record for record in SeqIO.parse(args.fasta, 'fasta')]
    # Calculate contig lengths
    contig_lengths = {record.id: len(record) for record in sequences}
    # Load the hallmark gene counts file
    hallmark_counts = pd.read_csv(args.hallmark_gene_summary, sep='\t')
    # Load the NCLDV score file
    NCLDV_assessment = pd.read_csv(args.ncldv_assessment, sep='\t')
    
    group1_reqs = check_requirements(group1, contig_lengths, hallmark_counts, single_copy_genes, NCLDV_assessment)
    group2_reqs = check_requirements(group2, contig_lengths, hallmark_counts, single_copy_genes, NCLDV_assessment)
    # Generate the summary DataFrame
    summary_df = create_summary_df(group1_reqs, group2_reqs, args.bin_id)

    # Save the summary DataFrame to a CSV file
    summary_df.to_csv(args.summary_output, index=False, sep='\t')

    processed_sequence_list = process_bin(summary_df, sequences, group1, group2)

    if len(processed_sequence_list) == 1:
        new_seq = processed_sequence_list[0]
        SeqIO.write(new_seq, f"{args.output}/{args.bin_id}.delineaged.fasta", "fasta")
    else:
        for i, new_seq in enumerate(processed_sequence_list):
            SeqIO.write(new_seq, f"{args.output}/{args.bin_id}_{i+1}.delineaged.fasta", "fasta")

if __name__ == "__main__":
    main()
