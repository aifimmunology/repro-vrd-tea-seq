from abc import ABC, abstractmethod
import polars as pl
from base64 import b64encode
from math import log10

# Example dataframe ---------------------------------------------------------------------------------

# shape: (46, 5)
# ┌───────────┬──────────┬──────────────┬─────────────────────────────────┬─────────────┐
# │ NES       ┆ adjP     ┆ -log10(adjP) ┆ Pathway                         ┆ Category    │
# │ ---       ┆ ---      ┆ ---          ┆ ---                             ┆ ---         │
# │ f64       ┆ f64      ┆ f64          ┆ str                             ┆ str         │
# ╞═══════════╪══════════╪══════════════╪═════════════════════════════════╪═════════════╡
# │ 1.810808  ┆ 0.370221 ┆ 0.431539     ┆ HALLMARK_CHOLESTEROL_HOMEOSTAS… ┆ not signif. │
# │ 1.698396  ┆ 0.370221 ┆ 0.431539     ┆ HALLMARK_APICAL_JUNCTION        ┆ not signif. │
# │ 1.455633  ┆ 0.370221 ┆ 0.431539     ┆ HALLMARK_HEME_METABOLISM        ┆ not signif. │
# │ 1.449694  ┆ 0.370221 ┆ 0.431539     ┆ HALLMARK_TGF_BETA_SIGNALING     ┆ not signif. │
# │ 1.423811  ┆ 0.370221 ┆ 0.431539     ┆ HALLMARK_PI3K_AKT_MTOR_SIGNALI… ┆ not signif. │
# │ …         ┆ …        ┆ …            ┆ …                               ┆ …           │
# │ -1.292415 ┆ 0.370221 ┆ 0.431539     ┆ HALLMARK_UV_RESPONSE_UP         ┆ not signif. │
# │ -1.309141 ┆ 0.407999 ┆ 0.389341     ┆ HALLMARK_MYC_TARGETS_V2         ┆ not signif. │
# │ -1.326192 ┆ 0.370221 ┆ 0.431539     ┆ HALLMARK_ALLOGRAFT_REJECTION    ┆ not signif. │
# │ -1.344837 ┆ 0.370221 ┆ 0.431539     ┆ HALLMARK_DNA_REPAIR             ┆ not signif. │
# │ -1.366121 ┆ 0.370221 ┆ 0.431539     ┆ HALLMARK_ESTROGEN_RESPONSE_EAR… ┆ not signif. │
# └───────────┴──────────┴──────────────┴─────────────────────────────────┴─────────────┘

# Another example
# shape: (890, 5)
# ┌───────────┬──────────┬──────────────┬───────────────┬─────────────┐
# │ NES       ┆ adjP     ┆ -log10(adjP) ┆ Pathway       ┆ Category    │
# │ ---       ┆ ---      ┆ ---          ┆ ---           ┆ ---         │
# │ f64       ┆ f64      ┆ f64          ┆ str           ┆ str         │
# ╞═══════════╪══════════╪══════════════╪═══════════════╪═════════════╡
# │ 2.34743   ┆ 0.001772 ┆ 2.751647     ┆ R-HSA-5663213 ┆ up          │
# │ 2.317871  ┆ 0.001772 ┆ 2.751647     ┆ R-HSA-6802946 ┆ up          │
# │ 2.317871  ┆ 0.001772 ┆ 2.751647     ┆ R-HSA-6802949 ┆ up          │
# │ 2.317871  ┆ 0.001772 ┆ 2.751647     ┆ R-HSA-6802955 ┆ up          │
# │ 2.317871  ┆ 0.001772 ┆ 2.751647     ┆ R-HSA-9649948 ┆ up          │
# │ …         ┆ …        ┆ …            ┆ …             ┆ …           │
# │ -1.651389 ┆ 0.349493 ┆ 0.456562     ┆ R-HSA-9734009 ┆ not signif. │
# │ -1.680046 ┆ 0.310211 ┆ 0.508342     ┆ R-HSA-1296071 ┆ not signif. │
# │ -1.686601 ┆ 0.335067 ┆ 0.474869     ┆ R-HSA-2172127 ┆ not signif. │
# │ -1.740465 ┆ 0.060577 ┆ 1.217695     ┆ R-HSA-909733  ┆ not signif. │
# │ -1.832037 ┆ 0.060577 ┆ 1.217695     ┆ R-HSA-2559584 ┆ not signif. │
# └───────────┴──────────┴──────────────┴───────────────┴─────────────┘
# NOTE: Notice how there is no pathway source column.

class GseaLoader(ABC):
    def __init__(
        self,
        p_col: str,
        es_col: str,
        store: str
    ):
        self.p_col = p_col
        self.es_col = es_col
        self.nlp_col = f"-log10({p_col})"

        # default cutoff values
        self.p_cutoff = -1*log10(0.05)
        self.es_cutoff = 0.1

        self.store = store

    @abstractmethod
    def load(self, df, labels, values, meta):
        pass

class GseaLoaderParquet(GseaLoader):

    def load(self, df, labels, values, meta):

        # import GSEA results
        df = pl.scan_parquet(df).select(labels + [self.es_col, self.p_col, "Pathway"])

        # filter results for covariate values and labels
        for label, value in zip(labels, values):
            df = df.filter(pl.col(label) == value)

        # assume pval needs negative log10 normalization
        df = df.with_columns(
            (-1 * pl.Expr.log10(pl.col(self.p_col) )).alias( self.nlp_col )
        )

        # label genes as up or down regulated
        df = df.with_columns(
            pl.when(pl.col(self.nlp_col).gt(self.p_cutoff))
            .then(
                pl.when(pl.col(self.es_col).lt(-1 * self.es_cutoff))
                    .then(pl.lit("down-regulated"))
                .when(pl.col(self.es_col).gt(self.es_cutoff))
                    .then(pl.lit("up-regulated"))
                .otherwise(pl.lit("not significant"))
            )
            .otherwise(pl.lit("not significant"))
            .alias("Category")
        )

        # collect dataframe without redundant information
        df = df.select([self.es_col, self.p_col, self.nlp_col, "Pathway", "Category"])

        # join dataframe with meta file to get Pathway Label
        meta_df = pl.scan_parquet(meta).select(["Pathway", "Pathway Label"])
        df = df.join(meta_df, on='Pathway', how='inner')

        # cast float64 to float32 then actualize DF
        df = df.cast({pl.Float64: pl.Float32})

        # export data
        serialized_df = b64encode(df.serialize()).decode('utf-8')

        return serialized_df
