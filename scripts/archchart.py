#!/usr/bin/env python3
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# The entry point function that based on final, filtered HMMsearch results
# generates HTML-formated depictions of all domain architectures found
# in HMMsearch results, also in regard to groups of domains. Arguments:
# --style   : CSS template for a class describing a domain group
# --tmpl    : a HTML template for the final output
# --colors  : a TSV file with data describing gradient star and end colors
#             for each group of domains
# --clstlen : a TSV file with cluster lenghts
# --hmmres  : pre-processed HMMsearch results in TSV format
# --output  : output FASTA file for extracted sequences
# USAGE:
# ./archchart.py --style CSS_TMPL   --tmpl HTML_TMPL    --colors DOM_COLORS \
#                --clstlen CLST_LEN --hmmres DOM_HMMOUT --output HTML_OUT

#-------------------------------------------------------------------------------
import argparse, sys
import pandas as pd
from os import linesep as eol

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --style   : CSS template for a class describing a domain group
       --tmpl    : a HTML template for the final output
       --colors  : a TSV file with data describing gradient star and end colors
                    for each group of domains
       --clstlen : a TSV file with cluster lenghts
       --hmmres  : pre-processed HMMsearch results in TSV format
       --output  : output FASTA file for extracted sequences
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--style', type=str, required=True,
        help='CSS template for a class describing a domain group')
    parser.add_argument('--tmpl', type=str, required=True,
        help='A HTML template for the final output')
    parser.add_argument('--colors', type=str, required=True,
        help='A TSV file with data describing gradient star and end colors ' +
             'for each group of domains')
    parser.add_argument('--clstlen', type=str, required=True,
        help='A TSV file with cluster lenghts')
    parser.add_argument('--hmmres', type=str, required=True,
        help='Final, filtered HMMsearch results in TSV format')
    parser.add_argument('--output', type=str, required=True,
        help='HTML output file with architecture charts')

    args = parser.parse_args()
    return args

def render_tab(df):
    '''A helper function that renders a HTML-formated table based on an input
       DataFrame. Arguments:
       df : an input DataFrame
       Returns:
       tab : string with generated HTML
    '''
    tab = f'<table><tr><th>{df.index.name}</th>'
    for name in df.columns:
        tab += f'<th>{name}</th>'
    tab += '</tr>' + eol
    for name, values in df.iterrows():
        tab += f'<tr><td>{name}</td>'
        for val in values:
            tab += f'<td class="val">{val:,}</td>'
        tab += '</tr>' + eol
    tab += '</table></br></br>'
    return tab

def main():
    '''The entry point function that based on final, filtered HMMsearch results
       generates HTML-formated depictions of all domain architectures found
       in HMMsearch results, also in regard to groups of domains.
    '''
    # Parse args and load domain group CSS style and final HTML templates.
    args = parse_args()
    with open(args.style) as f:
        style = f.read()
    with open(args.tmpl) as f:
        tmpl  = f.read()
    
    # Load tables assigning colors to group domains, add to it predefined colors
    # for domains unassigned to any group, next load preprocessed HMMsearch results.
    colors_df = pd.read_csv(args.colors, sep='\t')
    colors_df.loc[colors_df.shape[0]] = 'Unassigned', '#dddddd', '#cccccc'
    hmm_df    = pd.read_csv(args.hmmres, sep='\t')
    
    # Fill up CSS template with domain group names (class names) and colors
    # assigned to those domains in colors_df.
    styles = ''
    cols = 'group color1 color2'.split()
    for _, (group, color1, color2) in colors_df[cols].iterrows():
        styles += style.format(group=group, color1=color1, color2=color2) + eol
    # Create group_html and qname_html columns with domain names and accessions
    # flanked by <span> tags described by classes corresponding to domain groups.
    for group in colors_df['group']:
        for col in 'group', 'qname':
            hmm_df.loc[hmm_df['group'] == group, f'{col}_html'] = \
                f'<span class="{group}">' + \
                hmm_df.loc[hmm_df['group'] == group, col] + \
                '</span>'
    
    # In the preprocessed data rows have been already sorted by assembly accession,
    # target protein sequecne id/name and the star position of matched domain to
    # avoid assigning to a protein wrong architecture (domain order).
    # Having that done, group rows in respect to protein sequecne id/name,
    # join values of remaining colums (domain group and Pfam accession version)
    # with single dash. That creates for every target protein sequence a string
    # describing its domain architecture in regard to domain groups as well as
    # domain Pfam accession versions.
    # Replace single dashes with long dash character &#8212; for columns with
    # HTML-formated values.
    cols = 'tname clustid group qacc qname group_html qname_html'.split()
    agg_df = hmm_df[cols].groupby('clustid').agg(lambda vals: '--'.join(vals))
    for col in 'group_html', 'qname_html':
        agg_df.loc[:, col] = agg_df[col].str.replace('--', '&#8212;')
    
    # Read data on initial cluster lengths. Merge them with groupped HMMsearch
    # results on columns containing protein sequence id: clustid for clusters
    # and tname (target name) for HMMsearch results.
    print(agg_df.columns)
    clust_df = pd.read_csv(args.clstlen, sep='\t')
    print(clust_df.columns)
    merged_df = agg_df.merge(clust_df, how='left', left_on='clustid',
                             right_on='clustid')
    
    # Gather into a HTML table depictions of domain architectures and numbers of
    # linked to them representative as well as all sequences. Do that by Using
    # HTML-formatted values describing domain architectures (columns group_html
    # and qname_html), the number of sequences linked to them (non-redundant
    # cluster representatives, the number of sequences from HMMsearch results)
    # as well as the sum of lengths of clusters related to representative
    # sequecnes (all sequences linked to an architecture).
    content = ''
    for col, label in ('group_html', 'domain group'), ('qname_html', 'domain'):
        counts_df = merged_df[f'{col} length'.split()].groupby(col).sum()
        counts_df.sort_values('length', ascending=False, inplace=True)
        nr_counts = merged_df[f'{col} length'.split()].groupby(col).count()
        counts_df['NR Counts'] = nr_counts
        counts_df.rename({'length' : 'Counts'}, axis=1, inplace=True)
        counts_df.index.name = f'Architecture (by {label})'
        content += render_tab(counts_df) + eol
    # Put filled CSS template and rendered table into main HTML template and save.
    with open(args.output, 'w') as f:
        f.write(tmpl.format(styles=styles, content=content))

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

