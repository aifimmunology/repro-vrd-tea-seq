from io import BytesIO
from abc import ABC, abstractmethod
from base64 import b64encode, b64decode
from typing import Literal
from math import log10

import polars as pl

from dash import Input, Output, callback, no_update

from ..utils import serialize_lazyframe

# example inputs --------------------------------------------------------------

# example outputs -------------------------------------------------------------

# shape: (7_902, 7)
# ┌────────┬──────────┬───────────┬──────────────┬────────────┬────────────┬─────────────────┐
# │ Gene   ┆ padj     ┆ log2fc    ┆ -log10(padj) ┆ Foreground ┆ Background ┆ Category        │
# │ ---    ┆ ---      ┆ ---       ┆ ---          ┆ ---        ┆ ---        ┆ ---             │
# │ str    ┆ f32      ┆ f32       ┆ f32          ┆ str        ┆ str        ┆ str             │
# ╞════════╪══════════╪═══════════╪══════════════╪════════════╪════════════╪═════════════════╡
# │ A2M    ┆ 0.999989 ┆ -0.790189 ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# │ AAAS   ┆ 0.999989 ┆ 0.607844  ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# │ AACS   ┆ 0.999989 ┆ 0.45088   ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# │ AAGAB  ┆ 0.999989 ┆ 0.720984  ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# │ AAK1   ┆ 0.999989 ┆ 0.178421  ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# │ …      ┆ …        ┆ …         ┆ …            ┆ …          ┆ …          ┆ …               │
# │ ZXDC   ┆ 0.999989 ┆ 0.361337  ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# │ ZYG11B ┆ 0.999989 ┆ 1.38177   ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# │ ZYX    ┆ 0.999989 ┆ 0.267364  ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# │ ZZEF1  ┆ 0.999989 ┆ -0.126334 ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# │ ZZZ3   ┆ 0.999989 ┆ 0.009387  ┆ 0.000005     ┆ IL-18      ┆ PBS        ┆ not significant │
# └────────┴──────────┴───────────┴──────────────┴────────────┴────────────┴─────────────────┘

class DegByCovariatesReader(ABC):

    """
    Each derivative concrete class must have one .load() method.
    This class is never called in the main program because it doesn't
    have the right code to read the input files.
    """

    # document important attributes here
    # @property

    # NOTE: @abstractmethod means downstream concrete classes
    # NOTE: must provide implementation of the .load() method.
    @abstractmethod
    def load(self) -> pl.LazyFrame:
        """
        The polymorphic method .load() must contain one or 
        more dash callback and return a polars LazyFrame
        """
        pass

class DegByCovariatesFactory:

    """
    This class parses input files and file type to select the right class to read the file content.
    """

    # NOTE: @staticmethods are functions that don't rely on self, but can only invoked within this class.
    # NOTE: this helps code organization so we don't place this in utils.py that gets lost in the sauce.
    @staticmethod
    def select_reader(
        input_files:    str | list[str] | dict[str, str],
        file_type:      Literal["parquet", "csv"],
        **kwargs
    ):

        # check and instantiate parquet reader
        if isinstance(input_files, list) or isinstance(input_files, dict):
            if file_type == "parquet":

                # init concrete class
                deg_reader_object = DegByCovariatesReaderParquet()

                # invoke .load() method to register callbacks
                deg_reader_object.load(
                    files=input_files,
                    p_col=kwargs.get('p_col'),
                    es_col=kwargs.get('es_col'),
                    store=kwargs.get('store')
                )

                return None

        else:
            raise ValueError(f"No DESeq loader available for the following file: {input_files}")

class DegByCovariatesReaderParquet(DegByCovariatesReader):

    def load(
        self,
        files: dict,
        p_col: str,
        es_col: str,
        store: str
    ) -> pl.LazyFrame:

        # init vars
        self.p_col = p_col
        self.es_col = es_col
        self.nlp_col = f"-log10({p_col})"
        self.store = store
        t_adjp = pl.scan_parquet(files.get("t_adjp"))
        t_log2fc = pl.scan_parquet(files.get("t_log2fc"))

        # build the df
        @callback(
            Output(store, "data"),
            Input("gene-selection-store", "data"),
            Input("covariate-labels", "data"),
            Input("covariate-values", "data"),
            Input("dropdown_pval", "value"),
            Input("dropdown_es", "value"),
            Input("barcode", "data"),
            prevent_initial_call=True
        )
        def deg_by_covariates_pq(gene, labels, values, p_cutoff, es_cutoff, barcode):

            if barcode is None:
                return no_update

            if gene is None:
                gene="PDE4B"

            df = pl.concat([
                t_adjp.select([gene, "Barcode"]).rename({gene: self.p_col}),
                t_log2fc.select([gene, "Barcode"]).rename({gene: self.es_col}),
            ], how='align').with_columns(
                pl.col("Barcode").str.split_exact("_", len(labels)).struct.rename_fields(labels)
            ).unnest(
                "Barcode"
            ).filter(
                pl.col(labels[-1]).eq(values[-1])
            )

            # normalize pvalue to -log10 scale
            df = df.with_columns(
                (pl.Expr.log10(pl.col(self.p_col)).mul(-1)).alias(self.nlp_col)
            )

            df = df.with_columns(
                pl.concat_str(labels,  separator=": ").alias("Label")
            ).drop(
                labels
            )

            # insert significance cutoffs
            df = df.with_columns(
                pl.when(pl.col(self.nlp_col).gt(-1*log10(p_cutoff)))
                .then(
                    pl.when(pl.col(self.es_col).lt(-1 * es_cutoff))
                        .then(pl.lit("down-regulated"))
                    .when(pl.col(self.es_col).gt(es_cutoff))
                        .then(pl.lit("up-regulated"))
                    .otherwise(pl.lit("not significant"))
                )
                .otherwise(pl.lit("not significant"))
                .alias("Category")
            )

            # cast float64 to float32
            df = df.cast({pl.Float64: pl.Float32})

            # serialize data to pass onto next callback
            serialized_df = serialize_lazyframe(df)
            return serialized_df
