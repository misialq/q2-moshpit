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
        print(f'Fetching iPATH3 diagram for ids: {ids}')

    params = {
        'default_opacity': 0.5,
        'export_type': 'svg',
        'selection': '\n'.join([f'{x} {y} W15 1.0' for x, y in zip(ids, colors)])
    }

    response = requests.get(url=IPATH3_URL, params=params)

    with open(img_output_path, 'wb') as img:
        img.write(response.content)


def _process_slice(table: pd.DataFrame, slice_name: str) -> pd.DataFrame:
        """Process a data slice file into a DataFrame"""
        df = pd.read_csv(os.path.join(str(table), f'{slice_name}_slice.csv'), index_col='id')
        df = df[[col for col in df.columns if '(Intercept)' not in col]]
        df.columns = [slice_name]
        return df


def _generate_cmap(table: pd.DataFrame, output_dir: str) -> pd.Series:
    cmap = plt.get_cmap('coolwarm')
    min_lfc = table['lfc'].min()
    max_lfc = table['lfc'].max()

    def lfc_to_hex(x):
        r, g, b, _ = cmap((x - min_lfc) / (max_lfc - min_lfc))
        return '#%02x%02x%02x' % (int(r * 255), int(g * 255), int(b * 255))
    
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)

    norm = plt.Normalize(min_lfc, max_lfc)
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

    lfc_slice = _process_slice(differentials, 'lfc')
    q_slice = _process_slice(differentials, 'q_val')
    differentials = pd.concat([lfc_slice, q_slice], axis=1)

    # filter table by significance and LFC
    q_condition = (differentials['q_val'] < significance_threshold)
    lfc_condition = (differentials['lfc'] > log2_fold_change_threshold)
    differentials = differentials[q_condition & lfc_condition]
    
    # generate the color map
    differentials["color"] = _generate_cmap(differentials, output_dir)

    img_output_path = f'{output_dir}/pathway_map.svg'
    _fetch_ipath(differentials.index.tolist(), differentials['color'], img_output_path, verbose=True)

    TEMPLATES = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "assets",
        "exclusive"
    )
    shutil.copy(os.path.join(TEMPLATES, "index.html"), output_dir)
    shutil.copytree(os.path.join(TEMPLATES, "css"), os.path.join(output_dir, "css"))
    shutil.copytree(os.path.join(TEMPLATES, "js"), os.path.join(output_dir, "js"))

    templates = [os.path.join(TEMPLATES, "index.html"),]

    q2templates.render(templates, output_dir, context={"filenames": json.dumps(['pathway_map.svg'])})
