from abc import ABC
import polars as pl
from base64 import b64encode

from dash import Input, Output, callback, no_update

# example output ------------------------------------------------

# shape: (24, 5)
# ┌───────┬─────────┬──────────┬───────────────────────┬────────────┐
# │ Gene  ┆ Donor   ┆ Cytokine ┆ Cell Type             ┆ Expression │
# │ ---   ┆ ---     ┆ ---      ┆ ---                   ┆ ---        │
# │ str   ┆ str     ┆ str      ┆ str                   ┆ f64        │
# ╞═══════╪═════════╪══════════╪═══════════════════════╪════════════╡
# │ PDE4B ┆ Donor10 ┆ 4-1BBL   ┆ B-Intermediate/Memory ┆ 1.015893   │
# │ PDE4B ┆ Donor11 ┆ 4-1BBL   ┆ B-Intermediate/Memory ┆ 0.788048   │
# │ PDE4B ┆ Donor12 ┆ 4-1BBL   ┆ B-Intermediate/Memory ┆ 0.764743   │
# │ PDE4B ┆ Donor1  ┆ 4-1BBL   ┆ B-Intermediate/Memory ┆ 1.241975   │
# │ PDE4B ┆ Donor2  ┆ 4-1BBL   ┆ B-Intermediate/Memory ┆ 1.285135   │
# │ …     ┆ …       ┆ …        ┆ …                     ┆ …          │
# │ PDE4B ┆ Donor5  ┆ PBS      ┆ B-Intermediate/Memory ┆ 1.090341   │
# │ PDE4B ┆ Donor6  ┆ PBS      ┆ B-Intermediate/Memory ┆ 1.275891   │
# │ PDE4B ┆ Donor7  ┆ PBS      ┆ B-Intermediate/Memory ┆ 1.1734     │
# │ PDE4B ┆ Donor8  ┆ PBS      ┆ B-Intermediate/Memory ┆ 1.267148   │
# │ PDE4B ┆ Donor9  ┆ PBS      ┆ B-Intermediate/Memory ┆ 1.463739   │
# └───────┴─────────┴──────────┴───────────────────────┴────────────┘
# NOTE: the order of dataframe columns do matter! It should be in the following order:
# * Gene
# * Bulk label (Donor, Sample, Cluster, etc.)
# * Covariate labels
# * Expression

# abstract class (interface) ----------------------------------------------------------------------

# used to store common variables (arguments passed through concrete class)
class PseudoBulkLoader(ABC):
    def __init__(
        self,
        file_type:          str,
        meta_file:          str,
        norm_file:          str,
        feature_name:       str,

        perturbation:       bool,
        control:            str,
        bulk_label:         str,

        store:              str
    ):
        self.file_type      = file_type
        self.meta           = meta_file
        self.norm           = norm_file
        self.feature_name   = feature_name

        self.perturbation   = perturbation
        self.control        = control if perturbation else None
        self.bulk_label     = bulk_label
        self.store          = store

# concrete classes --------------------------------------------------------------------------------

# these classes inherit attributes from previous interface

class PseudoBulkLoaderParquet(PseudoBulkLoader):
    def __init__(self, file_type, meta_file, norm_file, feature_name, perturbation, control, bulk_label, store):
        super().__init__(file_type, meta_file, norm_file, feature_name, perturbation, control, bulk_label, store)

        # retrieve and save pb dataframe with input from gene click
        @callback(
            Output(self.store, "data"),
            Input("gene-selection-store", "data"),
            Input("barcode", "data"),
            Input("covariate-labels", "data"),
            Input("covariate-values", "data"),
            prevent_initial_call=True
        )
        def PseudoBulkLoaderParquet01(
            selected_feature,
            barcode,
            labels,
            values,
        ):

            if barcode is None:
                return no_update

            if selected_feature is None:
                selected_feature = "PDE4B"

            # remove any covariates with "Censor" in its labeling.
            pb_censor=False
            if pb_censor is False:
                labels = [ i for i in labels if "Censor" not in i ]
                values = [ i for i in values if "Censor" not in i ]

            # scan data
            pb_meta = pl.scan_parquet(self.meta)
            pb_norm = pl.scan_parquet(self.norm)

            # use the following chunk if perturbation experiment and control is provided
            if self.perturbation:
                if self.control:
                    control_condition="PBS"
                    for label, value in zip(labels, values):
                        pb_meta = pb_meta.filter(pl.col(label).is_in([value, control_condition]))

            barcodes = pb_meta.select("barcodes").collect().to_series().to_list()

            # re-shape data from wide to long format
            all_labels=[self.bulk_label] + labels

            # filter by selected gene
            df = pb_norm.filter(
                pl.col(self.feature_name) == selected_feature
            ).select(
                [self.feature_name] + barcodes
            ).unpivot(
                index=self.feature_name, variable_name="Covariates", value_name="Expression"
            ).with_columns(
                pl.col("Covariates").str.split_exact("_", len(labels)).struct.rename_fields(all_labels)
            ).unnest(
                "Covariates"
            )

            # cast float64 to float32
            df = df.cast({pl.Float64: pl.Float32})

            # serialize data to pass onto next callback
            serialized_df = b64encode(df.serialize()).decode('utf-8')
            return serialized_df


class GSVAScoreLoaderParquet(PseudoBulkLoader):
    def __init__(self, file_type, meta_file, norm_file, feature_name, perturbation, control, bulk_label, store):
        super().__init__(file_type, meta_file, norm_file, feature_name, perturbation, control, bulk_label, store)

        # retrieve and save pb dataframe with input from gene click
        @callback(
            Output(self.store, "data"),
            Input("dropdown_pathway_source", "value"),
            Input("pathway-selection-store", "data"),
            Input("barcode", "data"),
            Input("covariate-labels", "data"),
            Input("covariate-values", "data"),
            prevent_initial_call=True
        )
        def PseudoBulkLoaderParquet01(
            selected_source,
            selected_pathway,
            barcode,
            labels,
            values,
        ):

            if barcode is None:
                return no_update

            if selected_source is None:
                selected_source = "Hallmark"

            if selected_pathway is None:
                selected_pathway = "Inflammatory Response"

            cov = dict(zip(labels, values))

            # scan data
            pb_meta = pl.scan_parquet(self.meta)
            pb_norm = pl.scan_parquet(self.norm)

            pb_meta = pb_meta.filter(
                pl.col('Visit').is_in(
                    [cov['Foreground'], cov['Background']]
                ),
                pl.col('Cell Type') == cov['Cell Type']
            )

            barcodes = pb_meta.select("barcodes").collect().to_series().to_list()

            # re-shape data from wide to long format
            all_labels=[self.bulk_label] + labels

            source_pathway = selected_source + '_' + selected_pathway
            source_pathway = source_pathway.lower()

            # filter by selected gene
            df = pb_norm.with_columns(
                pl.col('Pathway').str.to_lowercase()
            ).filter(
                pl.col(self.feature_name) == source_pathway
            ).select(
                [self.feature_name] + barcodes
            ).unpivot(
                index = self.feature_name, 
                variable_name = 'barcodes',
                value_name = 'Score'
            ).join(
                pb_meta,
                how = 'left',
                on = 'barcodes'
            ).with_columns(
                pl.lit(selected_pathway).alias(self.feature_name)
            )

            # cast float64 to float32
            df = df.cast({pl.Float64: pl.Float32})

            # serialize data to pass onto next callback
            serialized_df = b64encode(df.serialize()).decode('utf-8')
            return serialized_df