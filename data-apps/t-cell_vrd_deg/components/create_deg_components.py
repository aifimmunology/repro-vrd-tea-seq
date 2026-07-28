from dash import dcc, html

import dash_ag_grid as dag

import polars as pl

import dash_bootstrap_components as dbc

from dash_iconify import DashIconify

import allen_dash_modules as adm

def create_deg_components(config: adm.ConfigSingleton) -> dict:

    def make_minimal_components(config: adm.ConfigSingleton) -> dict:

        # all finished components will be saved here and passed onto the DEG meta-class
        minimal_components = dict()

        # 1. build pvalue and enrichment score dropdowns
        dropdown_pval = create_pval_dropdown(id="dropdown_pval")
        dropdown_es = create_enrichment_score_dropdown(id="dropdown_es")
        minimal_components["dropdown_plot_filters"] = tuple([dropdown_pval, dropdown_es])

        # 2. build dropdowns for covariates
        dropdown_covariates = create_covariate_dropdowns(
            meta_df=config.data.deg.files.meta,
            covariates=config.experiment.covariates
        )
        minimal_components["dropdown_covariates"] = dropdown_covariates

        # 3. build dropdown for gene selection
        dropdown_gene = create_gene_dropdown(file=config.data.deg.files.adjp)
        minimal_components[dropdown_gene.id] = dropdown_gene

        # 4. create callback to use covariate values to parse DEG meta.parquet to get a barcode
        barcode_alert = dbc.Alert(
            children="",
            id="barcode_alert",
            is_open=False,
            color='danger',
            fade=True,
            dismissable=False
        )
        minimal_components["barcode_alert"] = barcode_alert

        # 5. volcano plot for deg...
        volcano_deg_plot = adm.VolcanoPlot(
            id="volcano_deg",
            section_title="Differential Gene Expression",
            p_col=config.data.deg.specifications.p_col,
            es_col=config.data.deg.specifications.es_col
        )
        minimal_components[volcano_deg_plot.id] = volcano_deg_plot

        # 6. define how to load deseq2 data
        adm.Deseq2Factory().select_reader(
            input_files=config.dict["data"]["deg"]["files"],
            file_type=config.data.deg.specifications.file_type,

            p_col=config.data.deg.specifications.p_col,
            lfc_col=config.data.deg.specifications.es_col,
            store=config.data.deg.specifications.store
        )

        # 6. table for deg...
        volcano_deg_grid = dag.AgGrid(id="volcano_deg_grid")
        minimal_components[volcano_deg_grid.id] = volcano_deg_grid

        # 7. table for covariates
        covariates_df_grid = dag.AgGrid(id="covariates_df_grid")
        minimal_components[covariates_df_grid.id] = covariates_df_grid

        # 8. sticky pathway feature
        sticky_pathway_checkbox = [
            dbc.Checkbox(
                id="sticky-pathway-checkbox",
                label="Enable sticky pathway feature",
                value=True,
                input_style={'width': '22px', 'height': '22px'},
                label_style={'fontSize': '1.25rem', 'marginLeft': '0.25rem'},
                style={
                    'display': 'inline-block',
                }
            ),
            DashIconify(
                id="sticky-pathway-icon",
                icon="material-symbols:help-outline",
                height=22,
                width=22,
                style={
                    'display': 'inline-block',
                    'marginBottom': '0.50rem',
                    'marginLeft': '0.25rem'
                }
            ),
            dbc.Popover(
                "Continue highlighting a selected pathway even \
                after changing dropdowns of independent variables",
                target="sticky-pathway-icon",
                trigger="hover",
                body=True
            )
        ]
        minimal_components["sticky_pathway_checkbox"] = sticky_pathway_checkbox

        # 9. GSEA rank toggle
        gsea_rank_toggle = [
            html.H5([
                'DEG ranking for GSEA',
                DashIconify(
                    id="gsea-radio-icon",
                    icon="material-symbols:help-outline",
                    className='iconify-in-text'
                )
            ]),
            dbc.Popover(
                "How should differentially expressed gene results be ranked \
                prior to gene set enrichment analyses?",
                target="gsea-radio-icon",
                trigger="hover",
                body=True
            ),
            dbc.RadioItems(
                options=[
                    {"label": "Direction & Adjusted P-value", "value": 'pval'},
                    {"label": "Wald statistic", "value": 'wald'}
                ],
                value='pval',
                id="gsea-rank-radio",
            )
        ]
        minimal_components["gsea_rank_toggle"] = gsea_rank_toggle

        return minimal_components
    
    def make_optional_components(config: adm.ConfigSingleton) -> dict:
        """
        optional_components is a dictionary that looks like this:
        
        {
            "optional": [gene_panel class, volcano_gsea, ...]
        }

        The minimal components dict gets updated with this dict to
        get assembled into the all_components.

        """

        optional_components = {}

        # include a gene panel if specified
        try:
            if config.dashboard.gene_panel:
                gene_panel = adm.GenePanel(id="gene_panel", gene_info_id="gene_info_panel")
                optional_components[gene_panel.id] = gene_panel
        except AttributeError as e:
            print(e)
            pass

        # include an MA plot if specified
        try:
            if config.data.deg.specifications.ma_plot:

                # initialize plot
                ma_plot = adm.MAPlot(
                    id="ma_plot",
                    section_title="MA Plot",
                    means_col=config.data.deg.specifications.means_col,
                    es_col=config.data.deg.specifications.es_col
                )

                # add component and store to output dict
                optional_components["ma-df-store"] = dcc.Store("ma-df-store")
                optional_components[ma_plot.id] = ma_plot

        except AttributeError as e:
            print(e)
            pass

        # include pseudobulk plot if specified
        try:
            if config.data.pseudobulk:
                pb = adm.BoxPlot(
                    id="boxplot_pb",
                    title="Pseudobulk Gene Expression"
                )

                # establish data-loader
                if config.data.pseudobulk.specifications.file_type == 'parquet':
                    adm.PseudoBulkLoaderParquet(
                        file_type=config.data.pseudobulk.specifications.file_type,
                        meta_file=config.data.pseudobulk.files.meta,
                        norm_file=config.data.pseudobulk.files.norm,
                        perturbation=config.experiment.perturbation,
                        control=config.experiment.control,
                        bulk_label=config.data.pseudobulk.specifications.bulk_label,
                        store=config.data.pseudobulk.specifications.store
                    )

                optional_components["pb-df-store"] = dcc.Store("pb-df-store")
                optional_components["pb_df_grid"] = dag.AgGrid(id="pb_df_grid")
                optional_components[pb.id] = pb

        except AttributeError as e:
            print("Omitting pseudobulk module from sc-deg-explorer.")
            pass

        try:
            if config.data.gsea:

                # since there are different sources of gsea, the app needs to generate the
                # volcano plot interactively. The p_col and es_col will be dynamically change
                # based on the input JSON, but the plot should still be generated and placed
                # into the app layout. This is a good example where the data source can change,
                # but Dash and ADM can still accommodate for these changes. therefore, the p_col
                # and es_col values are initially currently empty. 

                volcano_gsea = adm.VolcanoPlot(
                    id="volcano_gsea",
                    section_title="GSEA Pathways",
                    p_col="",
                    es_col=""
                )

                # establish choices of GSEA databases by looking into config file
                gsea_db_options = [ i for i in config.dict["data"]["gsea"].keys() if i != "specifications"]

                # create dd to specify GSEA source
                dropdown_pathway_source = create_pathway_source_dropdown(
                    options=gsea_db_options,
                    id="dropdown_pathway_source"
                )

                # create dd to specify GSEA pathway selection
                dropdown_pathway = create_pathway_dropdown(id="dropdown_pathway")
                optional_components[dropdown_pathway.id] = dropdown_pathway

                # place all created components to output dict
                optional_components["gsea-df-store"] = dcc.Store("gsea-df-store")
                optional_components["gsea_df_grid"] = dag.AgGrid(id="gsea_df_grid")
                optional_components[dropdown_pathway_source.id] = dropdown_pathway_source

                optional_components[volcano_gsea.id] = volcano_gsea

        except AttributeError as e:
            print("Omitting gsea module from sc-deg-explorer.")
            pass

        try:
            if config.data.gsva:

                # GSVA scores and stats are based on the same sets of pathways as GSEA results.
                # Therefore, dropdowns/selection are driven by the GSEA section above.
                # We'll build both a volcano plot and a box plot - one for stats and one for scores.

                # Volcano plot for GSVA stats
                volcano_gsva = adm.VolcanoPlot(
                    id="volcano_gsva",
                    section_title="GSVA Pathways",
                    p_col="",
                    es_col=""
                )

                gsva_stat_store = config.data.gsva.specifications.store

                # place volcano/stats components to output dict
                # plot data
                optional_components[gsva_stat_store] = dcc.Store(gsva_stat_store)
                # table data
                optional_components["gsva-stat_df_grid"] = dag.AgGrid(id="gsva-stat_df_grid")
                # volcano plot object
                optional_components[volcano_gsva.id] = volcano_gsva

                # Box plot for GSVA scores
                box_gsva = adm.BoxPlot(
                    id="boxplot_gsva",
                    title="GSVA Scores"
                )

                gsva_score_store = config.data.gsva_scores.specifications.store

                # establish data-loaders for scores
                if config.data.gsva_scores.specifications.file_type == 'parquet':
                    adm.GSVAScoreLoaderParquet(
                        file_type=config.data.gsva_scores.specifications.file_type,
                        meta_file=config.data.gsva_scores.files.meta,
                        norm_file=config.data.gsva_scores.files.norm,
                        feature_name="Pathway",
                        perturbation=config.experiment.perturbation,
                        control=config.experiment.control,
                        bulk_label=config.data.gsva_scores.specifications.bulk_label,
                        store=gsva_score_store
                    )

                # Stores for boxplot/score data
                # score data
                optional_components[gsva_score_store] = dcc.Store(gsva_score_store)
                # table data
                optional_components["gsva-score_df_grid"] = dag.AgGrid(id="gsva-score_df_grid")
                # box plot object
                optional_components[box_gsva.id] = box_gsva

        except AttributeError as e:
            print("Omitting gsva score module from sc-deg-explorer.")
            pass

        try:
            if config.data.deg.specifications.deg_by_covariates:

                # visualize a gene's expression across all covariates
                deg_by_covariates_plot = adm.VolcanoPlot(
                    id='deg_by_covariates_plot',
                    section_title='Gene Expression Across All Conditions',
                    p_col=config.data.deg.specifications.p_col,
                    es_col=config.data.deg.specifications.es_col
                )

                # select data reader for this plot
                adm.DegByCovariatesFactory().select_reader(
                    input_files=config.dict["data"]["deg"]["files"],
                    file_type=config.data.deg.specifications.file_type,

                    p_col=config.data.deg.specifications.p_col,
                    es_col=config.data.deg.specifications.es_col,
                    store="deg-by-covariates-df-store"
                )

                # create the modal that houses the plot
                deg_by_covariates_modal = dbc.Modal(
                    [
                        dbc.ModalHeader("", close_button=True),
                        dbc.ModalBody(deg_by_covariates_plot.ui),
                        dbc.ModalFooter([
                            html.P("This plot shows a gene's expression across every combination of experimental conditions."),
                            html.P("Click a point to change a set of experimental covariates.")
                            ]
                        )
                    ],
                        id="deg_by_covariates_modal",
                        is_open=False,
                        centered=True,
                        size='lg',
                        fade=True
                )

                # create the button that pulls up this plot
                deg_by_covariates_button = dbc.Button(
                    children="View Gene Expression for All Conditions",
                    id="deg_by_covariates_button",
                    outline=True,
                    color="primary"
                )

                # include stores and components into optional_components dict
                optional_components["deg-by-covariates-df-store"] = dcc.Store("deg-by-covariates-df-store")

                optional_components[deg_by_covariates_plot.id] = deg_by_covariates_plot
                optional_components[deg_by_covariates_modal.id] = deg_by_covariates_modal
                optional_components[deg_by_covariates_button.id] = deg_by_covariates_button

        except KeyError as e:
            print("Omitting deg's by covariates module from sc-deg-explorer.")
            pass

        try:
            if config.data.deg.specifications.summaries_plot:

                deg_summaries_plot = adm.DegSummariesPlot(
                    id="deg_summaries_plot",
                    section_title="Number of Differential Genes by Cell Type",
                    const_covariate = ["Cytokine", "Censoring"],
                    vary_covariates = ["Cell Type"]
                )

                deg_summaries_modal = dbc.Modal(
                    [
                        dbc.ModalHeader("", close_button=True),
                        dbc.ModalBody(deg_summaries_plot.ui),
                        dbc.ModalFooter([                            
                            html.P("This plot shows gene expression across all cell types in the specified set of conditions above."),
                            html.P("Clicking a bar will change the set of covariates and affect downstream plots.")
                            ]
                        )
                    ],
                        id="deg_summaries_modal",
                        is_open=False,
                        centered=True,
                        size='lg'
                )

                deg_summaries_button = dbc.Button(
                    "View Differential Genes for All Cell Types",
                    id="deg_summaries_button",
                    outline=True,
                    color="primary",
                )

                optional_components["deg-summaries-store"] = dcc.Store('deg-summaries-store')

                optional_components[deg_summaries_plot.id] = deg_summaries_plot
                optional_components[deg_summaries_modal.id] = deg_summaries_modal
                optional_components[deg_summaries_button.id] = deg_summaries_button

        except KeyError as e:
            print("Omitting summaries plot module from sc-deg-explorer.")
            pass

        return optional_components

    # invoke all sub-functions
    minimal_components = make_minimal_components(config)
    optional_components = make_optional_components(config)

    # merge minimal and optional dicts together. this feature only available in Python 3.9+
    all_components = minimal_components | optional_components

    return all_components

# supporting functions ----------------------------------------------------------------------------

def create_pval_dropdown(
    id="dropdown_pval",
    label="P-Value Cutoff",
    options=[ x/100.0 for x in range(5, 25, 5) ],
    **kwargs
    ):

    # create dropdown pvalue
    dropdown_pval = adm.DropDown(
        id=id,
        label=label,

        options=options,
        value=options[0],
        multi=False,
        clearable=False,
        **kwargs
    )

    return dropdown_pval

def create_enrichment_score_dropdown(
        id="dropdown_es",
        label="Log2FC Cutoff",
        options=[ x/100.0 for x in range(0, 25, 5) ][::-1],
        **kwargs
    ):

    dropdown_es = adm.DropDown(
        id=id,
        label=label,

        options=options,
        value=options[-1],
        multi=False,
        clearable=False,
        **kwargs
    )

    return dropdown_es

def create_covariate_dropdowns(meta_df, covariates: list):

    dropdown_covariates = list()
    meta_df = pl.read_parquet(meta_df)

    for covariate in covariates:

        # if there is a space in covariate, replace space with double underscore (will be removed during data parsing steps)
        dd_id = f"dropdown_{covariate}" if ' ' not in covariate else f"dropdown_{covariate.replace(' ', '__')}"

        # the options for a dropdown is the sorted values in the metadata df column of that covariate.
        options = sorted(meta_df[covariate].unique().to_list())

        # create one dd per covariate
        dd = adm.DropDown(
            id=dd_id,
            label=covariate,
            options=options,
            value=options[0],
            multi=False,
            clearable=False
        )

        # dropdown_covariates list
        dropdown_covariates.append(dd)

        if covariate == "Comparison Type":
            dd.dropdown.value = "Cross-Sectional (NDMM vs. Healthy)"
        if covariate == "Foreground":
            dd.dropdown.value = "Pre-Treatment (PreTx)"
            # dd.dropdown.value = "Healthy Controls"
        if covariate == "Background":
            dd.dropdown.value = "Healthy Controls"
            # dd.dropdown.value = "Pre-Treatment (PreTx)"
        if covariate == "Cell Type":
            dd.dropdown.value = "CD14 Mono Core"

    return tuple(dropdown_covariates)

def create_gene_dropdown(file, id="dropdown_gene", label="Gene", gene_col="Gene", **kwargs):

    # build gene dropdown
    options=pl.scan_parquet(file).select(gene_col).sort(by=gene_col).collect().to_series().to_list()

    dropdown_gene = adm.DropDown(
        id=id,
        label=label,

        options=options,
        multi=False,
        clearable=True,
        placeholder='Select a gene',
        **kwargs
    )

    # NOTE: While this component affects the value of gene-selection-store,
    # the logic to write to the gene store is located in the GeneSelectionMediator class.

    return dropdown_gene

def create_pathway_source_dropdown(options: list, id="dropdown_pathway_source", label="Pathway Source"):

    dropdown_pathway_source = adm.DropDown(
        id=id,
        label=label,
        options=options,
        value='Hallmark',
        multi=False,
        clearable=False
    )

    return dropdown_pathway_source

def create_pathway_dropdown(id="dropdown_pathway", label="Pathway"):

    # initally has None value but gets immediately initialized once the GSEA source is determined
    dropdown_pathway = adm.DropDown(
        id=id,
        label=label,
        value='Inflammatory Response',
        multi=False,
        clearable=True,
        placeholder='Select a Pathway'
    )

    return dropdown_pathway
