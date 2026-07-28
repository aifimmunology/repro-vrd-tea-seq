from io import BytesIO
from abc import ABC, abstractmethod
from base64 import b64encode, b64decode
from typing import Literal
from math import log10

import polars as pl

from dash import Input, Output, callback, no_update

import allen_dash_modules as adm

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

class Deseq2Reader(ABC):

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
    def load(self) -> None:
        """
        The polymorphic method .load() must contain one or 
        more dash callback and return a polars LazyFrame
        """
        pass

class Deseq2Factory:

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
                deseq2_reader_object = Deseq2ReaderParquet()

                # invoke .load() method to register callbacks
                deseq2_reader_object.load(
                    files=input_files,
                    p_col=kwargs.get('p_col'),
                    lfc_col=kwargs.get('lfc_col'),
                    store=kwargs.get('store')
                )

                return None

        else:
            raise ValueError(f"No DESeq loader available for the following file: {input_files}")

class Deseq2ReaderParquet(Deseq2Reader):

    def load(
        self,
        files: dict,
        p_col: str,
        lfc_col: str,
        store: str
    ) -> None:

        # how do I not repeat this section across all file formats,
        # while allowing the Factory to select the right concrete class?
        # it makes sense to display these attributes in the abstract class,
        # but since the abstract class never gets called then these attributes
        # are never realized. It's okay to repeat these across class as long as
        # it doesn't get out of hand.

        # set variables here
        self.pre_store = f"{store}-1"
        self.store = store

        self.p_col = p_col
        self.nlp_col = f"-log10({p_col})"
        self.lfc_col = lfc_col

        @callback(
            Output(self.pre_store, 'data'),
            Input('barcode', 'data'),
            prevent_initial_call=True
        )
        def DeseqReaderParquet01(
            barcode
        ):

            if barcode is None:
                return no_update

            # select relevant columns
            df = pl.concat([
                pl.scan_parquet(files["adjp"]).select(["Gene", barcode]).rename({barcode: 'padj'}),
                pl.scan_parquet(files["log2fc"]).select(["Gene", barcode]).rename({barcode: 'log2fc'})
            ], how='align')

            # assume pval needs negative log10 normalization
            df = df.with_columns(
                (-1 * pl.Expr.log10(pl.col(self.p_col) )).alias(self.nlp_col)
            )

            # Remove null values. They appear often when gene analysis was censored.
            df = df.drop_nulls()

            # insert Foreground/Background columns
            meta_df = pl.read_parquet(files["meta"]).filter(pl.col("barcodes") == barcode).select(["Foreground", "Background"])
            foreground = meta_df["Foreground"][0]
            background = meta_df["Background"][0]

            df = df.with_columns(
                pl.lit(foreground).alias("Foreground"),
                pl.lit(background).alias("Background")
            )

            # serialize data to pass onto next callback
            serialized_df = b64encode(df.serialize()).decode('utf-8')
            return serialized_df

        @callback(
            Output(self.store, 'data'),
            Input(self.pre_store, 'data'),
            Input('p-cutoff-store', 'data'),
            Input('es-cutoff-store', 'data'),
            prevent_initial_call=True
        )
        def DeseqReaderParquet02(df, p_cutoff, es_cutoff):

            # import data
            df = pl.LazyFrame.deserialize(BytesIO(b64decode(df)))

            # label genes as up or down regulated
            df = df.with_columns(
                pl.when(pl.col(self.nlp_col).gt(-1*log10(p_cutoff)))
                .then(
                    pl.when(pl.col(self.lfc_col).le(-1 * es_cutoff))
                        .then(pl.lit("down-regulated"))
                    .when(pl.col(self.lfc_col).ge(es_cutoff))
                        .then(pl.lit("up-regulated"))
                    .otherwise(pl.lit("not significant"))
                )
                .otherwise(pl.lit("not significant"))
                .alias("Category")
            )

            # cast float64 to float32 then actualize DF
            df = df.cast({pl.Float64: pl.Float32})

            # export data
            serialized_df = b64encode(df.serialize()).decode('utf-8')
            return serialized_df
