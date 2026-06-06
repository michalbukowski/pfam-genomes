#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Generates an HTML file depicting all domain architectures found in the
# integrated HMMsearch and SignalP results, including counts of representative
# and total sequences per architecture. Colours domain groups using CSS classes
# derived from a gradient color table. Arguments:
# --style   : CSS template for a class describing a domain group
# --tmpl    : HTML template for the final output
# --colors  : TSV file with gradient start and end colors for each domain group
# --clstlen : TSV file with cluster lengths (clustid and sequence counts)
# --allres  : integrated HMMsearch and SignalP results TSV
# --output  : output HTML file with domain architecture charts
# USAGE:
# ./archchart.py --style CSS_TMPL   --tmpl HTML_TMPL    --colors DOM_COLORS \
#                --clstlen CLST_LEN --allres INTEG_RES  --output HTML_OUT

#-------------------------------------------------------------------------------
import argparse, sys
import pandas as pd
from os import linesep as eol

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --style   : CSS template for a class describing a domain group
       --tmpl    : HTML template for the final output
       --colors  : TSV file with gradient start and end colors for each domain group
       --clstlen : TSV file with cluster lengths (clustid and sequence counts)
       --allres  : integrated HMMsearch and SignalP results TSV
       --output  : output HTML file with domain architecture charts
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
    parser.add_argument('--allres', type=str, required=True,
        help='Final HMMER and SignalP integrated results')
    parser.add_argument('--output', type=str, required=True,
        help='HTML output file with architecture charts')

    args = parser.parse_args()
    return args

def render_tab(df, colors_df):
    '''A helper function that renders a HTML-formated table based on an input
       DataFrame. Arguments:
       df   : an input DataFrame
       Returns:
       tab : string with generated HTML
    '''
    tab = f'<table><tr><th>{df.index.name}</th>'
    for name in df.columns[-2:]:
        tab += f'<th>{name}</th>'
    tab += '</tr>' + eol
    for _, (arch_group, arch_qname, *vals) in df.iterrows():
        arch = '&#8212;'.join([
            f'<span class="{group}">{qname}</span>'
            for group, qname in zip(arch_group, arch_qname)
        ])
        tab += f'<tr><td>{arch}</td>'
        for val in vals:
            tab += f'<td class="val">{val:,}</td>'
        tab += '</tr>' + eol
    tab += '</table></br></br>'
    return tab

def main():
    '''Loads integrated results and cluster lengths, builds CSS styles from the
       color table, aggregates domain hits per cluster into architecture strings,
       counts representative and total sequences per architecture, and writes
       an HTML file with a formatted architecture chart table.
    '''
    # Parse args and load domain group CSS style and final HTML templates.
    args = parse_args()
    with open(args.style) as f:
        style = f.read()
    with open(args.tmpl) as f:
        tmpl  = f.read()
    
    # Load tables assigning colors to domain groups and the final integrated
    # HMMsearch and SignalP results.
    colors_df = pd.read_csv(args.colors, sep='\t')
    all_df    = pd.read_csv(args.allres, sep='\t')

    # Fill up CSS template with domain group names (class names) and colors
    # assigned to those domains in colors_df.
    styles = ''
    cols = 'group color1 color2'.split()
    for _, (group, color1, color2) in colors_df[cols].iterrows():
        styles += style.format(group=group, color1=color1, color2=color2) + eol
    # Create qname_html column with domain names flanked by <span> tags described
    # by classes corresponding to domain groups.
    
    # Sort the final DataFrame with respect to domains locations (order domains).
    all_df.sort_values('tname env_from env_to'.split(), inplace=True)
    
    # Having the DataFrame sorted, group rows in respect to protein sequecne id/name,
    # join values of remaining colums (domain group and Pfam accession version)
    # with single dash. That creates for every target protein sequence a string
    # describing its domain architecture in regard to domain groups as well as
    # domain Pfam accession versions. Replace single dashes with long dash
    # character &#8212; for columns with HTML-formated values.
    # cols = 'tname clustid group qacc qname qname_html'.split()
    cols = 'clustid group qname'.split()
    agg_df = all_df[cols].groupby('clustid').agg(tuple)
    
    # Read data on initial cluster lengths. Merge them with groupped HMMsearch
    # results on columns containing protein sequence id: clustid for clusters
    # and tname (target name) for HMMsearch results.
    clust_df = pd.read_csv(args.clstlen, index_col='clustid', sep='\t')
    print(clust_df)
    print(agg_df)
    merged_df = agg_df.merge(clust_df, how='left', left_index=True, right_index=True)
    print(merged_df)
    # Gather into a HTML table depictions of domain architectures and numbers of
    # linked to them representative as well as all sequences. Do that by Using
    # HTML-formatted values describing domain architectures (column qname_html),
    # the number of sequences linked to them (non-redundant cluster representatives,
    # the number of sequences from HMMsearch results) as well as the sum of
    # lengths of clusters related to representative sequecnes (all sequences
    # linked to an architecture).
    content = ''
    merged_gb = merged_df['group qname length'.split()].groupby('group qname'.split())
    counts_df = merged_gb.sum().astype(int)
    counts_df.sort_values('length', ascending=False, inplace=True)
    nr_counts = merged_gb.count()
    counts_df['NR Count'] = nr_counts
    counts_df.reset_index(inplace=True)
    with_sp = counts_df['group'].apply(lambda groups: groups[0] == 'SP')
    counts_df = pd.concat([
        counts_df[  with_sp ],
        counts_df[ ~with_sp ]
    ])
    #counts_df.drop('qname', axis=1, inplace=True)
    counts_df.rename({'length' : 'Count'}, axis=1, inplace=True)
    counts_df.index.name = f'Domain architecture'
    content += render_tab(counts_df, colors_df) + eol
    # Put filled CSS template and rendered table into main HTML template and save.
    with open(args.output, 'w') as f:
        f.write(tmpl.format(styles=styles, content=content))

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

