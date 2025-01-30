# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import shutil
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

import pandas as pd
import q2templates
import requests
from q2_composition import DataLoafPackageDirFmt, FrictionlessCSVFileFormat

IPATH3_URL = 'https://pathways.embl.de/mapping.cgi'

def _fetch_ipath(ids: list, colors: list, img_output_path: str, verbose: bool = False):
    """Fetches an enriched pathways map from iPATH3 for given IDs."""
    # remove colon from EC names
    if ':' in ids[0]:
        ids = [x.replace(':', '') for x in ids]

    if verbose:
        print(f'Fetching iPATH3 diagram for {len(ids)} ids: {ids}')

    params = {
        'default_opacity': 0.2,
        'default_width': 3.0,
        'default_radius': 7.0,
        'export_type': 'svg',
        'background_color': '#ffffff',
        'default_color': '#aaaaaa',
        'map': 'metabolic',
        'keep_colors': True,
        'selection': '\n'.join([f'{x} {y} W12 1.0' for x, y in zip(ids, colors)])
    }

    response = requests.get(url=IPATH3_URL, params=params)

    with open(img_output_path, 'wb') as img:
        img.write(response.content)
    
    return params


def _process_slice(table: pd.DataFrame, slice_name: str) -> pd.DataFrame:
        """Process a data slice file into a DataFrame"""
        df = pd.read_csv(os.path.join(str(table), f'{slice_name}_slice.csv'), index_col='id')
        df = df[[col for col in df.columns if '(Intercept)' not in col]]
        return df


def _generate_cmap(table: pd.DataFrame, output_dir: str) -> pd.Series:
    cmap = plt.get_cmap('coolwarm')
    min_lfc = table['lfc'].min()
    max_lfc = table['lfc'].max()

    # Adjust the normalization to center around zero
    norm = plt.Normalize(vmin=min(min_lfc, 0), vmax=max(max_lfc, 0))

    def lfc_to_hex(x):
        # Map the log2 fold change to the colormap
        r, g, b, _ = cmap(norm(x))
        return '#%02x%02x%02x' % (int(r * 255), int(g * 255), int(b * 255))
    
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(sm, cax=ax, orientation='horizontal', ticklocation='top')
    cbar.set_label('log2 fold change', fontsize=12)

    plt.savefig(os.path.join(output_dir, 'legend.png'), bbox_inches='tight')

    return table['lfc'].apply(lfc_to_hex)


def map_pathways(
        output_dir: str,
        differentials: DataLoafPackageDirFmt,
        log2_fold_change_threshold: float = 0,
        significance_threshold: float = 0.05,
    ):
    # os.makedirs(os.path.join(output_dir, 'imgs'))

    # Process all comparisons, excluding "(Intercept)" columns
    lfc_slice = _process_slice(differentials, 'lfc')
    q_slice = _process_slice(differentials, 'q_val')

    comparisons = lfc_slice.columns.tolist()
    all_dfs = {}
    svg_filenames = []
    requests = {}
    for comparison in comparisons:
        lfc_df = lfc_slice.loc[:, comparison]
        lfc_df.name = 'lfc'
        q_df = q_slice.loc[:, comparison]
        q_df.name = 'q_val'

        comparison_df = pd.concat([lfc_df, q_df], axis=1)

        # filter table by significance and LFC
        q_condition = (comparison_df['q_val'] < significance_threshold)
        lfc_condition = (comparison_df['lfc'] > log2_fold_change_threshold)
        comparison_df = comparison_df[q_condition & lfc_condition]

        # generate the color map
        comparison_df["color"] = _generate_cmap(comparison_df, output_dir)

        all_dfs[comparison] = comparison_df

        img_output_path = f'{output_dir}/pathway_map_{comparison}.svg'
        svg_filenames.append(f'pathway_map_{comparison}.svg')
        params = _fetch_ipath(
            comparison_df.index.tolist(), comparison_df['color'],
            img_output_path, verbose=True
        )
        requests[comparison] = params

    TEMPLATES = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "assets",
        "exclusive"
    )
    shutil.copy(os.path.join(TEMPLATES, "index.html"), output_dir)
    shutil.copytree(os.path.join(TEMPLATES, "css"), os.path.join(output_dir, "css"))
    shutil.copytree(os.path.join(TEMPLATES, "js"), os.path.join(output_dir, "js"))

    templates = [os.path.join(TEMPLATES, "index.html"),]
    context = {
        "filenames": json.dumps(svg_filenames),
        "comparisons": json.dumps(comparisons),
        "request_body": json.dumps(requests),
    }

    q2templates.render(templates, output_dir, context=context)
