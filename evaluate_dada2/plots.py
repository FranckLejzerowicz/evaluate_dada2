# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plot_regressions(plots_pd, value_name, empties, mock_sam, f, r, p, pdf):
    g = sns.lmplot(
        data=plots_pd, x='mock (% reads)', y=value_name,
        hue='variable', col='comparison', height=4, aspect=0.8,
        facet_kws={'sharex': False, 'sharey': False})
    g.set_titles('{col_name}')
    plt.suptitle(
        '[%s-%s] Relative abundances of mock "%s" vs samples features (open-ref clust=%s)\n'
        '(%s %s of metadata samples do not shared any mock feature)' % (
            f, r, mock_sam, p, round(empties, 2), "%"), fontsize=14)
    plt.subplots_adjust(top=0.8)
    pdf.savefig(bbox_inches='tight')
    plt.close()


def make_heatmap_outputs(meta, stats_pd, pdf):
    """Make heatmaps for the DADA2 stats"""
    for value in ['passed filter', 'merged', 'non-chimeric']:
        name = 'percentage of input %s' % value
        fig, axes = plt.subplots(1, 2, figsize=(11, 3))
        for control in [0, 1]:
            sams = list(meta[meta['is_control'] == control].sample_name)
            cur_stats_pd = stats_pd[stats_pd['sample-id'].isin(sams)]
            stats_mean = cur_stats_pd.pivot_table(
                index=['forward'], columns=['reverse'],
                values=[name], aggfunc=np.mean)
            stats_sd = cur_stats_pd.pivot_table(
                index=['forward'], columns=['reverse'],
                values=[name], aggfunc=np.std)
            stats_full = round(
                stats_mean, 2).astype(str) + "\n(±" + round(
                stats_sd, 2).astype(str) + ")"
            stats_mean.columns = stats_mean.columns.droplevel()
            g = sns.heatmap(
                stats_mean, cmap='RdBu',
                annot=stats_full.values,
                ax=axes[control], fmt='')
            g.set_title('is control==%s (n=%s)' % (control, len(sams)))
        plt.suptitle(name, fontsize=14, fontweight="bold")
        plt.subplots_adjust(top=0.82)
        pdf.savefig(bbox_inches='tight')
        plt.close()


def make_heatmap_blast_asv(blast_in, pdf):
    """Parse the BLAST results to get the numbers of proper hits"""
    nqueries_pv = blast_in.pivot_table(
        index=['forward'],
        columns=['reverse'],
        values=['nqueries'])
    nqueries_pv.columns = nqueries_pv.columns.droplevel()
    fig, ax = plt.subplots(figsize=(4, 3))
    g = sns.heatmap(nqueries_pv, cmap='RdBu', annot=True, ax=ax)
    plt.suptitle('Number of ASVs in the trimmed mock sample',
                 fontsize=14, fontweight="bold")
    plt.subplots_adjust(top=0.82)
    pdf.savefig(bbox_inches='tight')
    plt.close()


def get_txts():
    txts = {
        'misclass': (
            0.7, 0.8,
            'False positives (misclassified vs "%s" [%s level])',
            "Features that do not match any expected features at the deepest level of classification (e.g., species level),\n" \
            "and usually represent either sample contaminants or sub-optimal bioinformatics pipelines (e.g., the presence of chimeric\n" \
            "sequences or the use of an overconfident taxonomic classifier)."
        ),
        'underclass': (
            0.65, 0.75,
            'False positives (underclassified vs "%s" [%s level])',
            "Observed features that match expected features, but are not classified to the expected taxonomic depth\n" \
            "(e.g., they are only classified to genus level but that genus classification is correct); these are often valid\n" \
            "features (i.e., not contaminants) but are not classified to the desired level either because of technical limitations\n" \
            "(e.g., sequences too short), degraded sequence quality, or sub-optimal methods."
        ),
        'false_neg': (
            0.72, 0.82,
            'False negatives vs "%s" (%s level)',
            "False negatives are features that were expected to be observed, but were not;\n" \
            "these can be compared to the false-positives to get an idea of what features may have been mis-/underclassified."
        ),
        'Observed Taxa': ("Obs features count", None),
        'Observed / Expected Taxa': ("Ratio of Obs:Exp features", None),
        'TAR': ("Taxon accuracy rate (TAR)",
                "Number of true pos features divided\n" \
                "by the total number of obs features\n" \
                "(true pos / (true pos + false pos))"),
        'TDR': ("Taxon detection rate (TDR)",
                "Number of true pos features divided\n" \
                "by the total number of exp features\n" \
                "(true pos / (true pos + false neg))"),
        'Bray-Curtis': ("Exp vs. Obs Bray-Curtis\ndissimilarity scores", None),
        'Jaccard': ("Exp vs. Obs Jaccard\ndistances scores", None),
        'r-squared': ("Exp vs. Obs linear\nregression r-squared value", None),
    }
    return txts


def make_heatmap_classifs(outs, txts, pdf):
    gb = ['type', 'sam']
    for typ in ['misclass', 'underclass', 'false_neg']:
        top, Y, suptitle, text = txts[typ]
        for (level, sam), gb_pd in outs[typ].groupby(gb, group_keys=False):
            n_p = gb_pd['p'].nunique()
            fig, axes = plt.subplots(1, n_p, figsize=(n_p * 12, 6))
            for pdx, (p, p_pd) in enumerate(gb_pd.groupby('p')):
                taxa_len = p_pd.pivot_table(
                    index=['f'], columns=['r'], values=['Taxon'],
                    aggfunc=lambda x: len(set([i for i in x if i != 'None'])))
                if level == 'taxo':
                    func = lambda x: '\n'.join(
                        sorted(set([i.split(';')[-1].strip() for i in x])))
                else:
                    func = lambda x: '\n'.join(
                        [';'.join(i) if len(set(x)) < 10 else '>10 IDs'
                         for i in np.array_split(sorted(set(x)), 4)])
                taxa_names = p_pd.pivot_table(
                    index=['f'], columns=['r'], values=['Taxon'],
                    aggfunc=func)
                taxa_relab = p_pd.pivot_table(
                    index=['f'], columns=['r'], values=['mock'],
                    aggfunc=sum)
                taxa_relab.columns = taxa_relab.columns.droplevel()
                taxa_names = taxa_len.astype(str) + '\n' + taxa_names.astype(
                    str)
                g = sns.heatmap(
                    taxa_relab, cmap='RdBu', fmt='', ax=axes[pdx],
                    annot=taxa_names.astype(str).values,
                    annot_kws={"fontsize": 6})
                g.set_title('Mock BLASTdb perc_ident = %s' % p, fontsize=12)
            plt.suptitle(suptitle % (sam, level), fontsize=20,
                         fontweight="bold")
            plt.text(.5, Y, text, transform=fig.transFigure,
                     horizontalalignment='center', fontsize=16)
            plt.subplots_adjust(top=top)
            pdf.savefig(bbox_inches='tight')
            plt.close()


def make_heatmap_stats(outs, txts, pdf):
    gb = ['type', 'level', 'sam', "p"]
    for (typ, level, sam, p), gb_pd in outs['res'].groupby(gb,
                                                           group_keys=False):
        if typ == 'taxo' and level == 1:
            continue
        fig, axes = plt.subplots(2, 3, figsize=(14, 7))

        obs, obsexp = 'Observed Taxa', 'Observed / Expected Taxa'
        obs_pv = gb_pd.pivot_table(index=['f'], columns=['r'], values=[obs])
        obs_pv.columns = obs_pv.columns.droplevel()
        obsexp_pv = round(100 * (
            gb_pd.pivot_table(index=['f'], columns=['r'], values=[obsexp])), 2)
        obsexp_pv.columns = obsexp_pv.columns.droplevel()
        both_pv = obs_pv.astype(str) + '\n (' + obsexp_pv.astype(str) + '%)'
        g = sns.heatmap(obs_pv, cmap='RdBu', fmt='', ax=axes[0, 0],
                        annot=both_pv.values, annot_kws={"fontsize": 6})
        g.set_title(txts[obsexp][0], fontsize=8)
        plt.text(.215, 0.885, txts[obs][0], transform=fig.transFigure,
                 horizontalalignment='center', fontsize=10, fontweight="bold")

        tar_pv = round(
            gb_pd.pivot_table(index=['f'], columns=['r'], values=['TAR']), 4)
        tar_pv.columns = tar_pv.columns.droplevel()
        g = sns.heatmap(tar_pv, cmap='RdBu', fmt='', ax=axes[0, 1],
                        annot=True, annot_kws={"fontsize": 6})
        g.set_title(txts['TAR'][1], fontsize=8)
        plt.text(.5, 0.92, txts['TAR'][0], transform=fig.transFigure,
                 horizontalalignment='center', fontsize=10, fontweight="bold")

        tdr_pv = round(
            gb_pd.pivot_table(index=['f'], columns=['r'], values=['TDR']), 4)
        tdr_pv.columns = tdr_pv.columns.droplevel()
        g = sns.heatmap(tdr_pv, cmap='RdBu', fmt='', ax=axes[0, 2],
                        annot=True, annot_kws={"fontsize": 6})
        g.set_title(txts['TDR'][1], fontsize=8)
        plt.text(.765, 0.92, txts['TDR'][0], transform=fig.transFigure,
                 horizontalalignment='center', fontsize=10, fontweight="bold")

        bc_pv = round(gb_pd.pivot_table(index=['f'], columns=['r'],
                                        values=['Bray-Curtis']), 4)
        bc_pv.columns = bc_pv.columns.droplevel()
        g = sns.heatmap(bc_pv, cmap='RdBu', fmt='', ax=axes[1, 0],
                        annot=True, annot_kws={"fontsize": 6})
        g.set_title(txts['Bray-Curtis'][0], fontsize=10, fontweight="bold")

        jc_pv = round(
            gb_pd.pivot_table(index=['f'], columns=['r'], values=['Jaccard']),
            4)
        jc_pv.columns = jc_pv.columns.droplevel()
        g = sns.heatmap(jc_pv, cmap='RdBu', fmt='', ax=axes[1, 1],
                        annot=True, annot_kws={"fontsize": 6})
        g.set_title(txts['Jaccard'][0], fontsize=10, fontweight="bold")

        r_pv = round(
            gb_pd.pivot_table(index=['f'], columns=['r'], values=['r-squared']),
            5)
        s_pv = round(
            gb_pd.pivot_table(index=['f'], columns=['r'], values=['Slope']), 2)
        p_pv = round(
            gb_pd.pivot_table(index=['f'], columns=['r'], values=['P value']),
            5)
        r_pv.columns = r_pv.columns.droplevel()
        s_pv.columns = s_pv.columns.droplevel()
        p_pv.columns = p_pv.columns.droplevel()
        annot_pv = '[s=' + s_pv.astype(str) + '\np=' + p_pv.astype(str) + ']'
        g = sns.heatmap(r_pv, cmap='RdBu', fmt='', ax=axes[1, 2],
                        annot=annot_pv.values, annot_kws={"fontsize": 5})
        plt.text(.765, 0.44, txts['r-squared'][0], transform=fig.transFigure,
                 horizontalalignment='center', fontsize=10, fontweight="bold")
        g.set_title('Linear regression stats [s=slope; p=p-value]', fontsize=8)

        if typ == 'asv':
            plt.suptitle(
                'Feature evaluation: mock "%s" vs ref (p=%s) ["%s" level]' % (
                    sam, p, typ), fontsize=15, fontweight="bold")
        else:
            plt.suptitle(
                'Feature evaluation: mock "%s" vs ref (p=%s) [%s level %s]' % (
                    sam, p, typ, level), fontsize=15, fontweight="bold")
        plt.subplots_adjust(top=0.85, hspace=0.525)
        pdf.savefig(bbox_inches='tight')
        plt.close()

