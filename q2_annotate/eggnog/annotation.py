# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
from io import StringIO

import pandas as pd
import qiime2

from q2_types.genome_data import (
    OrthologAnnotationDirFmt, Orthologs, SeedOrthologDirFmt, OrthologFileFmt
)
from q2_types.genome_data._deferred_setup._transformers import _annotations_to_dataframe, _is_valid_uuid4
from q2_types.reference_db import EggnogRefDirFmt
from q2_types.sample_data import SampleData


def _annotate_seed_orthologs_runner(
        seed_ortholog, eggnog_db, sample_label, output_loc,
        db_in_memory, num_cpus
):
    # at this point instead of being able to specify the type of target
    # orthologs, we want to annotate _all_.

    cmds = ['emapper.py', '-m', 'no_search', '--annotate_hits_table',
            str(seed_ortholog), '--data_dir', str(eggnog_db),
            '-o', str(sample_label), '--output_dir', str(output_loc),
            '--cpu', str(num_cpus)]
    if db_in_memory:
        cmds.append('--dbmem')

    subprocess.run(cmds, check=True, cwd=str(output_loc))


def _eggnog_annotate(
        eggnog_hits: SeedOrthologDirFmt,
        eggnog_db: EggnogRefDirFmt,
        db_in_memory: bool = False,
        num_cpus: int = 1
) -> OrthologAnnotationDirFmt:

    eggnog_db_fp = eggnog_db.path

    result = OrthologAnnotationDirFmt()

    # run analysis
    for relpath, obj_path in eggnog_hits.seed_orthologs.iter_views(
            OrthologFileFmt):
        sample_label = str(relpath).rsplit(r'.', 2)[0]

        _annotate_seed_orthologs_runner(seed_ortholog=obj_path,
                                        eggnog_db=eggnog_db_fp,
                                        sample_label=sample_label,
                                        output_loc=result,
                                        db_in_memory=db_in_memory,
                                        num_cpus=num_cpus)

    return result


def map_eggnog(
        ctx,
        eggnog_hits,
        eggnog_db,
        db_in_memory=False,
        num_cpus=1,
        num_partitions=None
):
    _eggnog_annotate = ctx.get_action("annotate", "_eggnog_annotate")
    collate_annotations = ctx.get_action(
        "types", "collate_ortholog_annotations"
    )

    if eggnog_hits.type <= SampleData[Orthologs]:
        partition_method = ctx.get_action("types", "partition_orthologs")
    else:
        raise NotImplementedError()

    (partitioned_orthologs, ) = partition_method(eggnog_hits, num_partitions)

    annotations = []
    for orthologs in partitioned_orthologs.values():
        (annotation, ) = _eggnog_annotate(
            orthologs, eggnog_db, db_in_memory, num_cpus
        )
        annotations.append(annotation)

    (collated_annotations, ) = collate_annotations(annotations)
    return collated_annotations


def _extract_generic(
        data: pd.DataFrame, column: str, transform_func: callable
) -> pd.Series:
    """
    Extracts a series from the annotation DataFrame and applies a
    transformation function to the specified column (annotation).

    Args:
        data (pd.DataFrame): The input DataFrame.
        column (str): The annotation column in the DataFrame to extract.
        transform_func (callable): The transformation function to apply
            to the extracted column.

    Returns:
        pd.Series: The extracted series with the transformation applied.
    """
    data = data.set_index("seed_ortholog", inplace=False)[column]
    data = data.apply(transform_func if data.notna().any() else lambda x: x)
    data = data.stack().reset_index(level=1, drop=True)
    data = data.value_counts()
    for char in ("", "-"):
        if char in data.index:
            data = data.drop(char, axis=0, inplace=False)
    return data


# this dictionary contains all the supported annotation types
# each value represents a tuple of:
# 1. original annotation column name (as it appears in the annotation table)
# 2. lambda function which will process values of that column and expand them
#   into a new series of values
extraction_methods = {
    "cog": ("COG_category", lambda x: pd.Series(list(x))),
    "kegg_ko": (
        "KEGG_ko",
        lambda x: pd.Series([i[3:] for i in x.split(",")])
    ),
    "kegg_pathway": (
        "KEGG_Pathway",
        lambda x: pd.Series([i for i in x.split(",") if i.startswith("map")])
    ),
    "kegg_module": ("KEGG_Module", lambda x: pd.Series(x.split(","))),
    "kegg_reaction": ("KEGG_Reaction", lambda x: pd.Series(x.split(","))),
    "brite": ("BRITE", lambda x: pd.Series(x.split(","))),
    "caz": ("CAZy", lambda x: pd.Series(x.split(",")))
}


def _filter(
        data: pd.DataFrame, max_evalue: float, min_score: float
) -> pd.DataFrame:
    data = data[(data["evalue"] <= max_evalue) & (data["score"] >= min_score)]
    if len(data) == 0:
        raise ValueError(
            "E-value/score filtering resulted in an empty table - "
            "please adjust your thresholds and try again."
        )
    return data


def extract_annotations(
        ortholog_annotations: OrthologAnnotationDirFmt,
        annotation: str,
        max_evalue: float = 1.0,
        min_score: float = 0.0
) -> pd.DataFrame:
    extract_method = extraction_methods.get(annotation)
    if not extract_method:
        raise NotImplementedError(
            f"Annotation '{annotation}' not supported."
        )
    else:
        col, func = extract_method

    annotations = []
    for _id, fp in ortholog_annotations.annotation_dict().items():
        annot_df = pd.read_csv(
            fp, sep="\t", skiprows=4, index_col=0
        )  # skip the first 4 rows as they contain comments
        annot_df = annot_df.iloc[:-3, :]  # remove the last 3 comment rows
        annot_df = _filter(annot_df, max_evalue, min_score)
        annot_df = _extract_generic(annot_df, col, func)
        annot_df.name = _id
        annotations.append(annot_df)

    result = pd.concat(annotations, axis=1).fillna(0).T
    result.index.name = "id"
    return result


def filter_annotations(
        ortholog_annotations: OrthologAnnotationDirFmt,
        query: str,
        max_evalue: float = 1.0,
        min_score: float = 0.0
) -> OrthologAnnotationDirFmt:
    annotations = ortholog_annotations.annotation_dict()
    filtered_annotations = OrthologAnnotationDirFmt()
    parsed = {}
    for _id, path in annotations.items():
        comments_start, comments_end, lines = [], [], []
        start = True
        with open(path, 'r') as f:
            for line in f:
                if line.startswith("##"):
                    if start:
                        comments_start.append(line)
                    else:
                        comments_end.append(line)
                else:
                    lines.append(line)
                    start = False

        df = pd.read_csv(StringIO('\n'.join(lines)), sep='\t', index_col=0)
        if _is_valid_uuid4(_id):
            df['MAG'] = _id
        else:
            df['Sample'] = _id

        # to satisfy QIIME2's particular requirements
        df.reset_index(drop=False, inplace=True)
        df.index = df.index.astype(str)
        df.index.rename('id', inplace=True)

        # filter by query
        metadata = qiime2.Metadata(df)
        ids = metadata.get_ids(where=query)
        if ids:
            filtered = metadata.filter_ids(ids).to_dataframe()
        else:
            filtered = pd.DataFrame()

        # filter by scores
        filtered = _filter(filtered, max_evalue, min_score)

        parsed[_id] = {
            'df': df, 'comments_start': comments_start,
            'comments_end': comments_end, 'filtered': filtered
        }

    for _id, values in parsed.items():
        fp = os.path.join(str(filtered_annotations), f"{_id}.emapper.annotations")
        with open(fp, 'w') as f:
            for comment in values['comments_start']:
                f.write(comment)

        # clean up the df
        for col in ['MAG', 'Sample']:
            if col in values['filtered']:
                values['filtered'].drop(col, axis=1, inplace=True)

        values['filtered'].to_csv(fp, mode='a', index=False, sep='\t')

        with open(fp, 'a') as f:
            for comment in values['comments_end']:
                f.write(comment)

    return filtered_annotations
