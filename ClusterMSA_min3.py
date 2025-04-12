import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np
import sys, os
from polyleven import levenshtein
from sklearn.cluster import DBSCAN
from utils import *

def plot_landscape(x, y, df, query_, plot_type):
    import seaborn as sns
    import matplotlib.pyplot as plt

    plt.figure(figsize=(20, 20))
    tmp = df.loc[df.dbscan_label == -1]
    plt.scatter(tmp[x], tmp[y], color='lightgray', marker='x', label='unclustered')

    tmp = df.loc[df.dbscan_label > 20]
    plt.scatter(tmp[x], tmp[y], color='yellow', label='other clusters')

    tmp = df.loc[(df.dbscan_label > 0) & (df.dbscan_label <= 20)]
    markers = ['o'] * 10 + ['s', 'P', 'D', 'X', '^', 'v', '<', '>', '8', 'H']
    palette = sns.color_palette('tab20', 20)
    marker_dict = {label: marker for label, marker in zip(range(1, 21), markers)}
    color_dict = {label: color for label, color in zip(range(1, 21), palette)}
    for label in range(1, 21):
        subset = tmp[tmp.dbscan_label == label]
        sns.scatterplot(x=x, y=y, data=subset, ax=plt.gca(), marker=marker_dict[label], 
                        color=color_dict[label], s=200, linewidth=0, label=f'Cluster {label}')

    plt.scatter(query_[x], query_[y], color='red', marker='*', s=450, label='Ref Seq')
    plt.legend(bbox_to_anchor=(1, 1), frameon=False)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.tight_layout()
    plt.savefig(args.o + "/" + args.keyword + '_' + plot_type + '.pdf', bbox_inches='tight')

if __name__=='__main__':
    p = argparse.ArgumentParser(description="""
    Cluster sequences in a MSA using DBSCAN algorithm and write .a3m file for each cluster.
    Assumes first sequence in fasta is the query sequence.
    H Wayment-Steele, 2022
    """)

    p.add_argument("keyword", action="store", help="Keyword to call all generated MSAs.")
    p.add_argument("-i", action='store', help='fasta/a3m file of original alignment.')
    p.add_argument("-o", action="store", help='name of output directory to write MSAs to.')
    p.add_argument("--n_controls", action="store", default=10, type=int, help='Number of control msas to generate (Default 10)')
    p.add_argument("--verbose", action='store_true', help='Print cluster info as they are generated.')
    p.add_argument("--resample", action='store_true', help='If included, will resample the original MSA with replacement before writing.')
    p.add_argument("--gap_cutoff", action='store', type=float, default=0.25, help='Remove sequences with gaps representing more than this frac of seq.')
    p.add_argument("--min_eps", action='store', default=3, type=float, help='Min epsilon value to scan for DBSCAN (Default 5).')
    p.add_argument("--eps_step", action='store', default=1.5, type=float, help='Step for epsilon scan for DBSCAN (Default 1.5).')
    p.add_argument("--min_samples", action='store', default=3, type=int, help='Default min_samples for DBSCAN (Default 3, recommended no lower than that).')
    p.add_argument("--log_dir", action="store", default=None, help="Directory to save the log file. Defaults to output directory.")

    args = p.parse_args()

    log_dir = args.log_dir if args.log_dir else args.o
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(args.o, exist_ok=True)

    log_file = os.path.join(log_dir, f"{args.keyword}.log")
    f = open(log_file, 'w')

    IDs, seqs = load_fasta(args.i)
    seqs = [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]

    df = pd.DataFrame({'SequenceName': IDs, 'sequence': seqs})
    query_ = df.iloc[:1]
    df = df.iloc[1:]

    if args.resample:
        df = df.sample(frac=1)

    L = len(df.sequence.iloc[0])
    df['frac_gaps'] = [x.count('-') / L for x in df['sequence']]

    former_len = len(df)
    df = df.loc[df.frac_gaps < args.gap_cutoff]
    new_len = len(df)
    lprint(args.keyword, f)
    lprint("%d seqs removed for containing more than %d%% gaps, %d remaining" % (former_len - new_len, int(args.gap_cutoff * 100), new_len), f)
    ohe_seqs = encode_seqs(df.sequence.tolist(), max_len=L)

    min_eps_values = [args.min_eps, args.min_eps + 5 * args.eps_step, args.min_eps + 10 * args.eps_step]
    all_eps_clusters = []

    for min_eps in min_eps_values:
        max_eps = min_eps + 10 * args.eps_step
        lprint(f'Testing for min_eps={min_eps}, max_eps={max_eps}', f)

        eps_test_vals = np.arange(min_eps, max_eps + args.eps_step, args.eps_step)
        for eps in eps_test_vals:
            testset = encode_seqs(df.sample(frac=0.25).sequence.tolist(), max_len=L)
            clustering = DBSCAN(eps=eps, min_samples=args.min_samples).fit(testset)
            n_clust = len(set(clustering.labels_))
            n_not_clustered = len(clustering.labels_[np.where(clustering.labels_ == -1)])
            lprint('%.2f\t%d\t%d' % (eps, n_clust, n_not_clustered), f)
            all_eps_clusters.append((eps, n_clust))

    eps_to_select, highest_clusters = max(all_eps_clusters, key=lambda x: x[1])
    lprint(f'Selected eps={eps_to_select:.2f} with {highest_clusters} clusters', f)

    clustering = DBSCAN(eps=eps_to_select, min_samples=args.min_samples).fit(ohe_seqs)
    lprint('Selected eps=%.2f' % eps_to_select, f)
    lprint("%d total seqs" % len(df), f)
    df['dbscan_label'] = clustering.labels_

    clusters = [x for x in df.dbscan_label.unique() if x >= 0]
    unclustered = len(df.loc[df.dbscan_label == -1])

    lprint('%d clusters, %d of %d not clustered (%.2f)' % (len(clusters), unclustered, len(df), unclustered / len(df)), f)

    avg_dist_to_query = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L for x in df.loc[df.dbscan_label == -1]['sequence'].tolist()])
    lprint('avg identity to query of unclustered: %.2f' % avg_dist_to_query, f)

    avg_dist_to_query = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L for x in df.loc[df.dbscan_label != -1]['sequence'].tolist()])
    lprint('avg identity to query of clustered: %.2f' % avg_dist_to_query, f)

    cluster_metadata = []
    for clust in clusters:
        tmp = df.loc[df.dbscan_label == clust]
        cs = consensusVoting(tmp.sequence.tolist())
        avg_dist_to_cs = np.mean([1 - levenshtein(x, cs) / L for x in tmp.sequence.tolist()])
        avg_dist_to_query = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L for x in tmp.sequence.tolist()])

        if args.verbose:
            print('Cluster %d consensus seq, %d seqs:' % (clust, len(tmp)))
            print(cs)
            print('#########################################')
            for _, row in tmp.iterrows():
                print(row['SequenceName'], row['sequence'])
            print('#########################################')

        tmp = pd.concat([query_, tmp], axis=0)
        cluster_metadata.append({
            'cluster_ind': clust,
            'consensusSeq': cs,
            'avg_lev_dist': '%.3f' % avg_dist_to_cs,
            'avg_dist_to_query': '%.3f' % avg_dist_to_query,
            'size': len(tmp)
        })
        write_fasta(tmp.SequenceName.tolist(), tmp.sequence.tolist(), outfile=args.o + '/' + args.keyword + '_' + "%03d" % clust + '.a3m')

    outfile = args.o + "/" + args.keyword + '_clustering_assignments.tsv'
    lprint('wrote clustering data to %s' % outfile, f)
    df.to_csv(outfile, index=False, sep='\t')

    metad_outfile = args.o + "/" + args.keyword + '_cluster_metadata.tsv'
    lprint('wrote cluster metadata to %s' % metad_outfile, f)
    metad_df = pd.DataFrame.from_records(cluster_metadata)
    metad_df.to_csv(metad_outfile, index=False, sep='\t')

    print(f'Saved this output to {log_file}')
    f.close()
