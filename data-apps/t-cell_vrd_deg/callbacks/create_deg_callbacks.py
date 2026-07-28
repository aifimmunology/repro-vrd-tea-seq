import sys
from io import BytesIO
from json import loads, dumps
from base64 import b64encode, b64decode

from dash import dcc, html, Input, Output, State, callback, ALL, no_update, ctx

import plotly.graph_objects as go
from plotly.express import scatter

from dash_iconify import DashIconify

import polars as pl

import allen_dash_modules as adm

from allen_dash_modules.utils import highlight_gene_members
from allen_dash_modules.utils import highlight_point_with_label
from allen_dash_modules.utils import highlight_selected_gene
from allen_dash_modules.utils import make_ag_grid_column_filters
from allen_dash_modules.utils import hex_to_rgba
from allen_dash_modules.utils import serialize_lazyframe
from allen_dash_modules.utils import deserialize_lazyframe
from allen_dash_modules.utils import patch_figure_theme

def create_deg_callbacks(config: adm.ConfigSingleton, components: dict):
    """
    The following sub-functions will define the 
    front and back-end behaviors of the app.
    """

    def create_layout_callbacks():
        """
        Required callbacks for the app layout to work correctly.
        """

        close_sidebar_icon = DashIconify(icon="bxs:left-arrow")
        open_sidebar_icon = DashIconify(icon="bxs:right-arrow")

        @callback(
            Output("sidebar", "className"),
            Output("main-content", "className"),
            Output("pull-tab", "className"),
            Output("pull-tab", "children"),
            Output("pull-tab-popover", "children"),
            Input("pull-tab", "n_clicks"),
            State("sidebar", "className"),
            prevent_initial_call=True
        )
        def toggle_sidebar(n, sidebar_class):
            """
            Open and close the sidebar.
            """
            if "collapsed" in sidebar_class:
                return (
                    "sidebar",
                    "main-content",
                    "pull-tab",
                    close_sidebar_icon,
                    "Collapse side panel"
                )
            else:
                return (
                    "sidebar collapsed",
                    "main-content expanded",
                    "pull-tab collapsed",
                    open_sidebar_icon,
                    "Expand side panel"
                )

        @callback(
            Output("barcode_alert", "is_open"),
            Output("barcode_alert", "children"),
            Input("barcode", "data")
        )
        def raise_barcode_error(barcode):

            if barcode is None:
                icon=DashIconify(icon="material-symbols:error-outline", height=26)
                header = html.Div(
                    [
                        html.Div(icon, style={'marginRight': '5px'}),
                        html.H4("Covariate Mismatch Error")
                    ],
                    style={'display': 'inline-flex'}
                )
                message_1 = html.P(
                    "Your experimental conditions did not lead to a valid analysis result! \
                    Please inspect the 'Covariates' table to find valid combinations of \
                    covariates."
                )
                message_2 = html.P(
                    "This error message will go away once a valid set of \
                    covariates is selected. Beware of interpreting any data \
                    interpretation at the moment!",
                    className='mb-0'
                )
                content=html.Div(
                    [
                        header,
                        html.Hr(),
                        message_1,
                        message_2
                    ]
                )

                return True, content
            else:
                return False, None
            

    def create_minimal_deg_callbacks(components: dict):

        # covariates --> covariate-values
        @callback(
            Output("covariate-values", "data"),
            [Input(component.id, "value") for component in components['dropdown_covariates']],
            prevent_initial_call='initial_duplicate'
        )
        def defineCovariateValues(*args):
            """
            Expected outputs are the values of the independent variables in tuple format.
            """
            return args

        # covariate-values --> barcode
        @callback(
            Output("barcode", "data"),
            Input("covariate-labels", "data"),
            Input("covariate-values", "data"),
            State("covariates-df-store", "data"),
            prevent_initial_call=True
        )
        def parse_barcode(labels, values, df):
            """
            Look through meta.parquet file to find a barcode.
            If successful, return a string; else return None.
            """

            # read meta.parquet file
            df = deserialize_lazyframe(df).collect()

            # filter dataframe
            if len(labels) == len(values):
                for lab, val in zip(labels, values):
                    df = df.filter(pl.col(lab) == val)
            else:
                sys.exit("Mismatch between number of covariates and values! " \
                "Check for errors in the pattern-matching callbacks")

            # count number of rows post filtering
            n_barcodes = df.shape[0]

            # return barcode if only one match
            if n_barcodes == 1:
                return "_".join(values)

            # return None if no match
            elif n_barcodes == 0:
                return None

            # if there are multiple matching barcodes, raise an error!
            elif n_barcodes >1:
                sys.exit("ERROR: multiple matching covariate barcodes detected! " \
                "Please inspect experimental setup and configuration.")
                return None

        # pval dropdown --> pval-cutoff-store
        @callback(
            Output("p-cutoff-store", "data"),
            Input("dropdown_pval", "value")
        )
        def return_pval_cutoff_store(value):
            """
            Expected output is the p-value cutoff as a float
            """
            return value

        # log2fc dropdown --> log2fc-cutoff-store
        @callback(
            Output("es-cutoff-store", "data"),
            Input("dropdown_es", "value")
        )
        def return_es_cutoff_store(value):
            """
            Expected output is the enrichment score cutoff as a float
            """
            return value

        # DESeq2 data-loading schema if a callback and not an object.

        # DESeq2 volcano callback to plot figure
        volcano_deg_plot = components["volcano_deg"]

        @callback(
            Output(volcano_deg_plot.id, 'figure'),
            Input('volcano-deg-df', 'data'),
            Input('toggle-dark-mode', 'value'),
            Input('pathway-selection-store', 'data'),
            Input('gene-selection-store', 'data'),
            State('pathway-gene-store', 'data'),
            prevent_initial_call=True
        )
        def plot_volcano_deg(
            df,
            light,
            pathway_selection_store,
            gene,
            pathway_gene_store
        ):

            # patch figure theme if light toggle
            if ctx.triggered_id == "toggle-dark-mode":
                fig = patch_figure_theme(light)
                return fig

            # import data
            df = deserialize_lazyframe(df).collect()

            # create plot
            fig = volcano_deg_plot.plot(
                data_frame = df,
                template="plotly_white" if light else "plotly_dark",
                custom_data=["Gene", volcano_deg_plot.p_col]
            )

            # format hover template
            fig.update_traces(
                hovertemplate=
                    '<b>%{customdata[0]}</b>' + '<br><br>' +
                    f'<b>{volcano_deg_plot.p_col}: </b>' + '%{customdata[1]:.2e}' + '<br>' +
                    f'<b>{volcano_deg_plot.nlp_col}: </b> ' + '%{y:.2f}' + '<br>' +
                    f'<b>{volcano_deg_plot.es_col}: </b> ' + '%{x:.2f}' + '<br>' +
                    "<extra></extra>",
            )

            # add gene regulation direction labels
            fig.add_annotation(text=f'Higher in {df["Foreground"][0]}', showarrow=False, x=1, y=1.05, xanchor='right', yanchor='top', xref='paper', yref='paper', font=dict(size=16))
            fig.add_annotation(text=f'Higher in {df["Background"][0]}', showarrow=False, x=0, y=1.05, xanchor='left', yanchor='top', xref='paper', yref='paper', font=dict(size=16))

            # highlight gene members if pathway was selected
            if pathway_selection_store is not None:

                # de-serialize encoded json
                pathway_gene_store = loads(b64decode(pathway_gene_store.encode('utf-8')).decode('utf-8'))

                gene_members = pathway_gene_store[pathway_selection_store]
                highlight_gene_members(
                    df = df,
                    fig = fig,
                    feature = "Gene",
                    genes = gene_members,
                    x_col=volcano_deg_plot.es_col,
                    y_col=volcano_deg_plot.nlp_col,
                    color_map=volcano_deg_plot.point_colors
                )

            # highlight gene with respect to gene-selection-store
            if gene is not None:
                highlight_selected_gene(
                    df = df,
                    fig = fig,
                    feature = "Gene",
                    gene = gene,
                    x_col=volcano_deg_plot.es_col,
                    y_col=volcano_deg_plot.nlp_col,
                    color_map=volcano_deg_plot.point_colors,
                    title=None
                )

            return fig

        # another callback to configure behavior of covariates_df_grid
        @callback(
            [Output(comp.id, "value") for comp in components["dropdown_covariates"]],
            Input("covariates_df_grid", "cellClicked"),
            Input("covariates-df-store", "data"),
            Input("covariate-labels", "data"),
            prevent_initial_call=True
        )
        def covariate_df_change_covariate_values(grid_click, df, labels):

            if grid_click:

                # get rowId from user click. not rowIndex because
                # that changes value if the table is sorted or filtered
                row_index = int(grid_click["rowId"])
                
                # import dataframe as lazyframe
                df = pl.LazyFrame.deserialize(BytesIO(b64decode(df)))

                # slice lazyframe by index, then get covariate values given labels
                df = df.with_row_index().filter(pl.col("index") == row_index).select(labels)

                # collect first (and only) row in the dataframe into a tuple.
                covariate_tuple = df.collect().row(0)

                return covariate_tuple

            else:
                return no_update

        # DEG AgGrid Table
        @callback(
            Output("volcano_deg_grid", "rowData"),
            Output("volcano_deg_grid", "columnDefs"),
            Output("volcano_deg_grid", "className"),
            Output("volcano_deg_grid", "getRowStyle"),
            Output("volcano_deg_grid", "scrollTo"),
            Input("volcano-deg-df", "data"),
            Input("gene-selection-store", "data"),
            Input("toggle-dark-mode", "value"),
            Input("tertiary-tabs", "active_tab"),
            Input("volcano_deg_grid", "columnState"),
            Input("volcano_deg_grid", "filterModel"),
            prevent_initial_call=True
        )
        def deg_df_grid_callback(
            df,
            gene,
            light,
            active_tab,
            column_state,
            filter_model
        ):

            # import data
            df = pl.LazyFrame.deserialize(BytesIO(b64decode(df)))

            # define all column fields
            fields = make_ag_grid_column_filters(df)

            # table theme
            className = "ag-theme-alpine" if light else "ag-theme-alpine-dark"

            if gene is None:
                scroll_index = None
                row_style = None
            else:

                # filter for feature
                filtered_df = df.with_row_index().filter(
                    pl.col("Gene") == gene
                ).collect()

                if filtered_df.is_empty():
                    scroll_index = None
                    row_style = None
                else:

                    # make scroll index
                    data_row = df.filter(
                        pl.col("Gene") == gene
                    )

                    scroll_index = {
                        'data': data_row.collect().to_dicts()[0],
                        'rowPosition': 'middle'
                    }

                    # make row style
                    row_style = {
                        "styleConditions": [{
                            "condition": f"params.data.Gene === '{gene}'",
                            "style": {"backgroundColor": "rgba(25, 118, 210, 0.40)"}
                        }]
                    }
            # trigger callback if active tab
            if active_tab == "DEG Table":
                return df.collect().to_dicts(), fields, className, row_style, scroll_index

            # return tbl properties
            return df.collect().to_dicts(), fields, className, row_style, scroll_index

        @callback(
            Output("covariates_df_grid", "rowData"),
            Output("covariates_df_grid", "columnDefs"),
            Output("covariates_df_grid", "className"),
            Output("covariates_df_grid", "getRowStyle"),
            Output("covariates_df_grid", "scrollTo"),
            Input("covariates-df-store", "data"),
            Input("barcode", "data"),
            Input("toggle-dark-mode", "value"),
            Input("tertiary-tabs", "active_tab"),
            Input("covariates_df_grid", "columnState"),
            Input("covariates_df_grid", "filterModel"),
            prevent_initial_call=True
        )
        def covariates_df_grid_callback(
            df,
            barcode,
            light,
            active_tab,
            column_state,
            filter_model
        ):

            # import data
            df = pl.LazyFrame.deserialize(BytesIO(b64decode(df)))

            # no update if barcode is not formed well
            if barcode is None:
                return no_update

            # table theme
            className = "ag-theme-alpine" if light else "ag-theme-alpine-dark"

            # make row style
            row_style = {
                "styleConditions": [{
                    "condition": f"params.data.barcodes === '{barcode}'",
                    "style": {"backgroundColor": "rgba(25, 118, 210, 0.40)"}
                }]
            }

            # make scroll index
            data_row = df.filter(
                pl.col("barcodes") == barcode
            )

            scroll_index = {
                'data': data_row.collect().to_dicts()[0],
                'rowPosition': 'middle'
            }

            # define all column fields
            fields = make_ag_grid_column_filters(df, exclude=["barcodes"])

            # trigger callback if active tab
            if active_tab == "Covariates":
                return df.collect().to_dicts(), fields, className, row_style, scroll_index

            # return tbl properties
            return df.collect().to_dicts(), fields, className, row_style, scroll_index

        # another callback to configure behavior of covariates_df_grid
        @callback(
            [Output(comp.id, "value", allow_duplicate=True) for comp in components["dropdown_covariates"]],
            Input("covariates_df_grid", "cellClicked"),
            Input("covariates-df-store", "data"),
            Input("covariate-labels", "data"),
            prevent_initial_call=True
        )
        def covariate_df_change_covariate_values(grid_click, df, labels):

            if grid_click:

                # get rowId from user click. not rowIndex because
                # that changes value if the table is sorted or filtered
                row_index = int(grid_click["rowId"])
                
                # import dataframe as lazyframe
                df = pl.LazyFrame.deserialize(BytesIO(b64decode(df)))

                # slice lazyframe by index, then get covariate values given labels
                df = df.with_row_index().filter(pl.col("index") == row_index).select(labels)

                # collect first (and only) row in the dataframe into a tuple.
                covariate_tuple = df.collect().row(0)

                return covariate_tuple

            else:
                return no_update

        # gene selection mediator
        gene_selection_mediator = adm.GeneSelectionMediator()
        gene_selection_mediator.mediate(
            deg_plot=volcano_deg_plot,
            dropdown_gene=components["dropdown_gene"]
        )

    def create_optional_deg_callbacks(config: adm.ConfigSingleton, components: dict):
        
        if components.get("ma_plot"):

            ma_plot=components["ma_plot"]

            # config is only necessary to get file locations.
            # should they be defined upon component creation?
            # aka added as a class attribute...

            @callback(
                Output("ma-df-store", "data"),
                Input("barcode", "data"),
                Input("es-cutoff-store", "data"),
                prevent_initial_call=True
            )
            def ma_df(barcode, es_cutoff):

                if barcode is None:
                    return no_update

                log2fc = pl.scan_parquet(config.data.deg.files.log2fc)
                means = pl.scan_parquet(config.data.deg.files.means)
                meta_df = pl.scan_parquet(config.data.deg.files.meta).filter(pl.col("barcodes")==barcode)
                foreground = meta_df.select("Foreground").first().collect().item()
                background = meta_df.select("Background").first().collect().item()

                # read gene means and log2fc
                df = pl.concat([
                    means.select(["Gene", barcode]).rename({barcode: ma_plot.means_col}),
                    log2fc.select(["Gene", barcode]).rename({barcode: ma_plot.es_col}),
                ], how='align').drop_nulls()
                
                # filter away highly-expressing genes
                df = df.filter(pl.col("means") < 4000)

                # insert foreground/background columns
                df = df.with_columns(
                    pl.lit(foreground).alias("Foreground"),
                    pl.lit(background).alias("Background"),
                )

                # cast float64 to float32
                df = df.cast({pl.Float64: pl.Float32})

                # insert significance cutoffs
                df = df.with_columns(
                    pl.when(pl.col(ma_plot.es_col).lt(-1 * es_cutoff))
                        .then(pl.lit("down-regulated"))
                    .when(pl.col(ma_plot.es_col).gt(es_cutoff))
                        .then(pl.lit("up-regulated"))
                    .otherwise(pl.lit("not significant"))
                    .alias("Category")
                )

                # serialize data to pass onto next callback
                serialized_df = serialize_lazyframe(df)
                return serialized_df

            @callback(
                Output(ma_plot.id, 'figure'),
                Input("ma-df-store", 'data'),
                Input("es-cutoff-store", "data"),
                Input("gene-selection-store", "data"),
                Input("pathway-gene-store", "data"),
                Input("pathway-selection-store", "data"),
                Input("toggle-dark-mode", "value"),
                prevent_initial_call=True
            )
            def make_ma_plot(df, es_cutoff, gene, pathway_gene_store, pathway, light):

                # patch figure theme if light toggle
                if ctx.triggered_id == "toggle-dark-mode":
                    fig = patch_figure_theme(light)
                    return fig

                # import data
                df = deserialize_lazyframe(df).collect()

                # create plot
                fig = ma_plot.plot(
                    data_frame = df,
                    es_cutoff=es_cutoff,
                    template="plotly_white" if light else "plotly_dark",
                    custom_data=["Gene", "means", "log2fc"]
                )

                fig.add_annotation(text=f'Higher in {df["Foreground"][0]}', showarrow=False, x=1, y=1.05, xanchor='right', yanchor='top', xref='paper', yref='paper', font=dict(size=16))
                fig.add_annotation(text=f'Higher in {df["Background"][0]}', showarrow=False, x=1, y=0, xanchor='right', yanchor='bottom', xref='paper', yref='paper', font=dict(size=16))

                if pathway is not None:

                    # de-serialize encoded json
                    pathway_gene_store = loads(b64decode(pathway_gene_store.encode('utf-8')).decode('utf-8'))

                    gene_members = pathway_gene_store[pathway]
                    highlight_gene_members(
                        df = df,
                        fig = fig,
                        feature = "Gene",
                        genes = gene_members,
                        x_col=ma_plot.means_col,
                        y_col=ma_plot.es_col,
                        color_map=ma_plot.point_colors
                    )

                if gene is not None:
                    highlight_selected_gene(
                        df = df,
                        fig = fig,
                        feature = "Gene",
                        gene = gene,
                        x_col="means",
                        y_col="log2fc",
                        color_map=ma_plot.point_colors,
                        title=None
                    )

                return fig

        if components.get("boxplot_pb"):

            boxplot_pb = components["boxplot_pb"]

            @callback(
                Output(boxplot_pb.id, "figure"),
                Input("pb-df-store", "data"),
                Input("gene-selection-store", "data"),
                Input("toggle-dark-mode", "value"),
                prevent_initial_call=True
            )
            def plot_pb(df, gene, light):

                # import data
                df = deserialize_lazyframe(df)

                # cast float64 to float32 then actualize DF
                df = df.cast({pl.Float64: pl.Float32}).collect()

                # create title
                gene=df['Gene'].unique()[0]
                plot_title=f"Selected Gene: <b>{gene}</b>"

                # define templates
                template="plotly_white" if light else "plotly_dark"

                # Group by donor and collect their cytokine values
                # TODO: optimize plotting strategy
                fig = go.Figure()

                fig = boxplot_pb.plot(
                    df=df,
                    title=plot_title,
                    template=template,
                    condition="Cytokine",
                    perturbation=True,
                    perturbation_subject_label="Donor",
                    perturbation_control_label=config.experiment.control
                )

                return fig

            @callback(
                Output("pb_df_grid", "rowData"),
                Output("pb_df_grid", "columnDefs"),
                Output("pb_df_grid", "className"),
                Output("pb_df_grid", "getRowStyle"),
                Output("pb_df_grid", "scrollTo"),
                Input("pb-df-store", "data"),
                Input("barcode", "data"),
                Input("toggle-dark-mode", "value"),
                Input("tertiary-tabs", "active_tab"),
                prevent_initial_call=True
            )
            def pb_df_grid_callback(
                df,
                barcode,
                light,
                active_tab
            ):

                # import data
                df = pl.LazyFrame.deserialize(BytesIO(b64decode(df)))

                # define all column fields
                fields = make_ag_grid_column_filters(df, exclude="barcodes")

                # no update if barcode is not formed well
                if barcode is None:
                    return no_update

                # table theme
                className = "ag-theme-alpine" if light else "ag-theme-alpine-dark"

                # make row style
                row_style = None

                # adjust tbl scroll
                scroll_index = None

                # trigger callback if active tab
                if active_tab == "Pseudo-Bulk Table":
                    return df.collect().to_dicts(), fields, className, row_style, scroll_index

                # return tbl properties
                return df.collect().to_dicts(), fields, className, row_style, scroll_index

        # establish callbacks for gsea pathway analysis
        if components.get("volcano_gsea"):

            volcano_gsea = components["volcano_gsea"]

            # allow source dropdown to change the options in the pathways
            @callback(
                Output("dropdown_pathway", "options"),
                Output("pathway-label-store", "data"),
                Output("pathway-gene-store", "data"),
                Input("dropdown_pathway_source", "value")
            )
            def update_pathway_options_given_pathway_source(source):
                """
                This callback accomplishes the following tasks:
                    1. Sets the GSEA database options via the config file.

                    2. Reads a meta.parquet file to establish
                        which gene members are in a pathway via dict, which looks like...
                        {
                            "HALLMARK_IL2_STAT5_SIGNALING": ('ABCB1', 'ADAM19', 'AGER', ...),
                            "HALLMARK_ADIPOGENESIS": ('ADD1', 'AIFM3', 'ANKH', ...)
                        }

                    3. Reads a a meta.parquet file to establish
                        relation between pathway code and pathway label, which looks like...
                        {
                            "HALLMARK_IL2_STAT5_SIGNALING": "IL-2/STAT5 Signaling",
                            "HALLMARK_ADIPOGENESIS": "Adipogenesis"
                        }

                These operations are placed inside a callback because GSEA
                sources can change during app operation, which differ from other
                data-loading schemes that are determined from app startup.
                """

                # this dict stores pathway and gene members (for #2 above)
                pathway_gene_dict = dict()

                # this dict stores pathway and label (for #3 above)
                pathway_label_dict = dict()

                # read the config file to pull up specifications and files
                gsea_files = config.dict["data"]["gsea"][source]["files"]
                gsea_meta = pl.scan_parquet(gsea_files["meta"])

                # use this loop to populate dicts for tasks #2 and #3
                for pathway, label, genes in gsea_meta.select(
                    ["Pathway", "Pathway Label", "Pathway Genes"]
                ).collect().rows():
                    if label is not None and genes is not None:
                        pathway_label_dict[pathway] = label
                        pathway_gene_dict[pathway] = tuple(genes.split(";"))

                # serialize dictionaries and send as output
                pathway_label_serialized = b64encode(dumps(pathway_label_dict).encode('utf-8')).decode('utf-8')
                pathway_gene_serialized = b64encode(dumps(pathway_gene_dict).encode('utf-8')).decode('utf-8')

                # options for gsea dropdown
                gsea_dropdown = sorted(tuple(pathway_label_dict.values()))

                return gsea_dropdown, pathway_label_serialized, pathway_gene_serialized

            @callback(
                Output("gsea-df-store", 'data'),
                Input("dropdown_pathway_source", "value"),
                Input("gsea-rank-radio", "value"),
                Input("covariate-labels", "data"),
                Input("covariate-values", "data"),
                Input("barcode", "data"),
                prevent_initial_call=True
            )
            def build_volcano_gsea_df_pq(
                source,
                rank,
                labels,
                values,
                barcode
            ):
                if barcode is None:
                    return no_update

                # NOTE: we cannot re-initialize a volcano_gsea because the previously made one
                # with an empty p_col and es_col variables is more 'globally' scoped than below.
                # We can only modify the attributes of the previously-made object and have the
                # changes propagate.

                # volcano_gsea = adm.VolcanoPlot(
                #     id="volcano_gsea",
                #     section_title="GSEA Pathways",
                #     p_col=config.dict["data"]["gsea"][source]["specifications"]["p_col"],
                #     es_col=config.dict["data"]["gsea"][source]["specifications"]["es_col"]
                # )

                # given a pre-made volcano plot, adjust its attributes.
                volcano_gsea.p_col = config.dict["data"]["gsea"][source]["specifications"]["p_col"]
                volcano_gsea.es_col = config.dict["data"]["gsea"][source]["specifications"]["es_col"]
                volcano_gsea.nlp_col = f"-log10({volcano_gsea.p_col})"

                # dynamically determine the class that will load the data.
                file_type=config.dict["data"]["gsea"][source]["specifications"]["file_type"]
                
                if rank == 'pval':
                    results_file=config.dict["data"]["gsea"][source]["files"]["log_results"]
                elif rank == 'wald':
                    results_file=config.dict["data"]["gsea"][source]["files"]["wald_results"]
                meta_file=config.dict["data"]["gsea"][source]["files"]["meta"]

                p_cutoff = config.dict["data"]["gsea"]["specifications"]["log10_p_cutoff"]
                es_cutoff = config.dict["data"]["gsea"]["specifications"]["es_cutoff"]

                if file_type == "parquet":

                    # define input columns to read
                    df_columns = labels + [volcano_gsea.es_col, volcano_gsea.p_col, "Pathway"]
                    if "Foreground" not in labels:
                        df_columns.append("Foreground")
                    if "Background" not in labels:
                        df_columns.append("Background")

                    # import gsea results
                    df = pl.scan_parquet(results_file).select(df_columns)
                    # define barcode column and filter for it
                    df = df.with_columns(
                        pl.concat_str(labels, separator="_").alias("barcode")
                    ).filter(
                        pl.col("barcode") == barcode
                    )

                    # assume pval needs negative log10 normalization
                    df = df.with_columns(
                        (-1 * pl.col(volcano_gsea.p_col).log(10)).alias( volcano_gsea.nlp_col )
                    )

                    # label genes as up or down regulated
                    df = df.with_columns(
                        pl.when(pl.col(volcano_gsea.nlp_col).gt(p_cutoff))
                        .then(
                            pl.when(pl.col(volcano_gsea.es_col).lt(-1 * es_cutoff))
                                .then(pl.lit("down-regulated"))
                            .when(pl.col(volcano_gsea.es_col).gt(es_cutoff))
                                .then(pl.lit("up-regulated"))
                            .otherwise(pl.lit("not significant"))
                        )
                        .otherwise(pl.lit("not significant"))
                        .alias("Category")
                    )

                    # collect dataframe with minimal info
                    df = df.select([volcano_gsea.es_col, volcano_gsea.p_col, volcano_gsea.nlp_col, "Pathway", "Category", "Foreground", "Background"])

                    # join df with meta file to get Pathway Label
                    meta_df = pl.scan_parquet(meta_file).select(["Pathway", "Pathway Label"])
                    df = df.join(meta_df, on='Pathway', how='inner')

                    # cast float64 to float32 then actualize DF
                    df = df.cast({pl.Float64: pl.Float32})

                    # export data
                    gsea_df = serialize_lazyframe(df)

                    # return the dataframe
                    return gsea_df

            @callback(
                Output(volcano_gsea.id, 'figure'),
                Input('gsea-df-store', 'data'),
                Input('toggle-dark-mode', 'value'),
                Input('pathway-selection-store', 'data'),
                State('pathway-label-store', 'data'),
                prevent_initial_call=True
            )
            def plot_volcano_gsea(
                df,
                light,
                pathway_selection_store,
                pathway_label_store
            ):

                # patch figure theme if light toggle
                if ctx.triggered_id == "toggle-dark-mode":
                    fig = patch_figure_theme(light)
                    return fig

                # import data
                df = deserialize_lazyframe(df).collect()

                # plot figure
                fig = volcano_gsea.plot(
                    df,
                    template="plotly_white" if light else "plotly_dark",
                    custom_data=["Pathway Label", "adjP"]
                )

                fig.update_traces(
                    hovertemplate=
                    '<b> ' + "%{customdata[0]}" + '</b><br>' +
                    "<br>" +
                    f"{volcano_gsea.p_col} = " + "%{customdata[1]:.2e}<br>" +
                    f"{volcano_gsea.nlp_col} = " + "%{y:.2f}<br>" +
                    f"{volcano_gsea.es_col} = " + "%{x:.3f}<br>" +
                    "<extra></extra>"
                )

                fig.add_annotation(text=f'Higher in {df["Foreground"][0]}', showarrow=False, x=1, y=1.05, xanchor='right', yanchor='top', xref='paper', yref='paper', font=dict(size=16))
                fig.add_annotation(text=f'Higher in {df["Background"][0]}', showarrow=False, x=0, y=1.05, xanchor='left', yanchor='top', xref='paper', yref='paper', font=dict(size=16))

                if pathway_selection_store is not None:

                    # de-serialize pathway label store
                    pathway_label_store = loads(b64decode(pathway_label_store.encode('utf-8')).decode('utf-8'))

                    highlight_point_with_label(
                        df = df,
                        fig = fig,
                        feature_col = "Pathway Label",
                        feature_value = pathway_label_store[pathway_selection_store],
                        x_col=volcano_gsea.es_col,
                        y_col=volcano_gsea.nlp_col,
                        color_map=volcano_gsea.point_colors,
                        title=True
                    )

                return fig

            # GSEA AgGrid Table
            @callback(
                Output("gsea_df_grid", "rowData"),
                Output("gsea_df_grid", "columnDefs"),
                Output("gsea_df_grid", "className"),
                Output("gsea_df_grid", "getRowStyle"),
                Output("gsea_df_grid", "scrollTo"),
                Input("gsea-df-store", "data"),
                Input("pathway-selection-store", "data"),
                Input("toggle-dark-mode", "value"),
                Input("tertiary-tabs", "active_tab"),
                Input("gsea_df_grid", "columnState"),
                Input("gsea_df_grid", "filterModel"),
                prevent_initial_call=True
            )
            def gsea_df_grid_callback(
                df,
                pathway,
                light,
                active_tab,
                column_state,
                filter_model
            ):

                # import data
                df = pl.LazyFrame.deserialize(BytesIO(b64decode(df)))

                # clean data
                df = df.sort(
                    "Pathway Label"
                ).filter(
                    pl.col("Pathway Label").is_not_null()
                )

                # table theme
                className = "ag-theme-alpine" if light else "ag-theme-alpine-dark"

                if pathway is None:
                    scroll_index = None
                    row_style = None
                else:
                    # filter for feature
                    filtered_df = df.with_row_index().filter(
                        pl.col("Pathway") == pathway
                    ).collect()

                    if filtered_df.is_empty():
                        scroll_index = None
                        row_style = None
                    else:

                        # make scroll index
                        data_row = df.filter(
                            pl.col("Pathway") == pathway
                        )

                        scroll_index = {
                            'data': data_row.collect().to_dicts()[0],
                            'rowPosition': 'middle'
                        }

                        # make row style
                        row_style = {
                            "styleConditions": [{
                                "condition": f"params.data.Pathway === '{pathway}'",
                                "style": {"backgroundColor": "rgba(25, 118, 210, 0.40)"}
                            }]
                        }

                # sort and define all column fields
                df = df.select(
                    ["Pathway", "Pathway Label", pl.all().exclude(["Pathway", "Pathway Label"])]
                )

                fields = make_ag_grid_column_filters(df, exclude=["Pathway"])

                # trigger callback if active tab
                if active_tab == "Pathway Table":
                    return df.collect().to_dicts(), fields, className, row_style, scroll_index

                # return tbl properties
                return df.collect().to_dicts(), fields, className, row_style, scroll_index
        
        # establish callbacks for GSVA pathway analysis
        # No pathway drop-down logic - that's all handled by the GSEA panels, above.
        if components.get("volcano_gsva"):

            volcano_gsva = components["volcano_gsva"]

            @callback(
                Output("gsva-stat-df-store", 'data'),
                Input("dropdown_pathway_source", "value"),
                Input("covariate-labels", "data"),
                Input("covariate-values", "data"),
                Input("barcode", "data"),
                prevent_initial_call=True
            )
            def build_volcano_gsva_df_pq(
                source,
                labels,
                values,
                barcode
            ):

                if barcode is None:
                    return no_update

                # NOTE: we cannot re-initialize a volcano_gsea because the previously made one
                # with an empty p_col and es_col variables is more 'globally' scoped than below.
                # We can only modify the attributes of the previously-made object and have the
                # changes propagate.

                # volcano_gsea = adm.VolcanoPlot(
                #     id="volcano_gsea",
                #     section_title="GSEA Pathways",
                #     p_col=config.dict["data"]["gsea"][source]["specifications"]["p_col"],
                #     es_col=config.dict["data"]["gsea"][source]["specifications"]["es_col"]
                # )

                # given a pre-made volcano plot, adjust its attributes.
                volcano_gsva.p_col = config.dict["data"]["gsva"][source]["specifications"]["p_col"]
                volcano_gsva.es_col = config.dict["data"]["gsva"][source]["specifications"]["es_col"]
                volcano_gsva.nlp_col = f"-log10({volcano_gsva.p_col})"

                # dynamically determine the class that will load the data.
                file_type=config.dict["data"]["gsva"][source]["specifications"]["file_type"]

                results_file=config.dict["data"]["gsva"][source]["files"]["results"]
                meta_file=config.dict["data"]["gsea"][source]["files"]["meta"]

                p_cutoff = config.dict["data"]["gsva"]["specifications"]["log10_p_cutoff"]
                es_cutoff = config.dict["data"]["gsva"]["specifications"]["es_cutoff"]

                if file_type == "parquet":

                    # define input columns to read
                    df_columns = labels + [volcano_gsva.es_col, volcano_gsva.p_col, "Pathway"]
                    if "Foreground" not in labels:
                        df_columns.append("Foreground")
                    if "Background" not in labels:
                        df_columns.append("Background")

                    # import gsea results
                    df = pl.scan_parquet(results_file).select(df_columns)

                    # define barcode column and filter for it
                    df = df.with_columns(
                        pl.concat_str(labels, separator="_").alias("barcode")
                    ).filter(
                        pl.col("barcode") == barcode
                    )

                    # assume pval needs negative log10 normalization
                    df = df.with_columns(
                        (-1 * pl.col(volcano_gsva.p_col).log(10)).alias( volcano_gsva.nlp_col )
                    )

                    # label genes as up or down regulated
                    df = df.with_columns(
                        pl.when(pl.col(volcano_gsva.nlp_col).gt(p_cutoff))
                        .then(
                            pl.when(pl.col(volcano_gsva.es_col).lt(-1 * es_cutoff))
                                .then(pl.lit("down-regulated"))
                            .when(pl.col(volcano_gsva.es_col).gt(es_cutoff))
                                .then(pl.lit("up-regulated"))
                            .otherwise(pl.lit("not significant"))
                        )
                        .otherwise(pl.lit("not significant"))
                        .alias("Category")
                    )

                    # collect dataframe with minimal info
                    df = df.select([volcano_gsva.es_col, volcano_gsva.p_col, volcano_gsva.nlp_col, "Pathway", "Category", "Foreground", "Background"])

                    # join df with meta file to get Pathway Label
                    meta_df = pl.scan_parquet(meta_file)#.select(["Pathway", "Pathway Label"])
                    
                    df = df.join(meta_df, on='Pathway', how='inner')

                    # cast float64 to float32 then actualize DF
                    df = df.cast({pl.Float64: pl.Float32})

                    # export data
                    gsva_df = serialize_lazyframe(df)

                    # return the dataframe
                    return gsva_df

            @callback(
                Output(volcano_gsva.id, 'figure'),
                Input('gsva-stat-df-store', 'data'),
                Input('toggle-dark-mode', 'value'),
                Input('pathway-selection-store', 'data'),
                State('pathway-label-store', 'data'),
                prevent_initial_call=True
            )
            def plot_volcano_gsva(
                df,
                light,
                pathway_selection_store,
                pathway_label_store
            ):

                # patch figure theme if light toggle
                if ctx.triggered_id == "toggle-dark-mode":
                    fig = patch_figure_theme(light)
                    return fig

                # import data
                df = deserialize_lazyframe(df).collect()

                # plot figure
                fig = volcano_gsva.plot(
                    df,
                    template="plotly_white" if light else "plotly_dark",
                    custom_data=["Pathway Label", "adjP"]
                )

                fig.update_traces(
                    hovertemplate=
                    '<b> ' + "%{customdata[0]}" + '</b><br>' +
                    "<br>" +
                    f"{volcano_gsva.p_col} = " + "%{customdata[1]:.2e}<br>" +
                    f"{volcano_gsva.nlp_col} = " + "%{y:.2f}<br>" +
                    f"{volcano_gsva.es_col} = " + "%{x:.3f}<br>" +
                    "<extra></extra>"
                )

                fig.add_annotation(text=f'Higher in {df["Foreground"][0]}', showarrow=False, x=1, y=1.05, xanchor='right', yanchor='top', xref='paper', yref='paper', font=dict(size=16))
                fig.add_annotation(text=f'Higher in {df["Background"][0]}', showarrow=False, x=0, y=1.05, xanchor='left', yanchor='top', xref='paper', yref='paper', font=dict(size=16))

                if pathway_selection_store is not None:

                    # de-serialize pathway label store
                    pathway_label_store = loads(b64decode(pathway_label_store.encode('utf-8')).decode('utf-8'))

                    highlight_point_with_label(
                        df = df,
                        fig = fig,
                        feature_col = "Pathway Label",
                        feature_value = pathway_label_store[pathway_selection_store],
                        x_col=volcano_gsva.es_col,
                        y_col=volcano_gsva.nlp_col,
                        color_map=volcano_gsva.point_colors,
                        title=True
                    )

                return fig

            # GSVA AgGrid Table
            @callback(
                Output("gsva-stat_df_grid", "rowData"),
                Output("gsva-stat_df_grid", "columnDefs"),
                Output("gsva-stat_df_grid", "className"),
                Output("gsva-stat_df_grid", "getRowStyle"),
                Output("gsva-stat_df_grid", "scrollTo"),
                Input("gsva-stat-df-store", "data"),
                Input("pathway-selection-store", "data"),
                Input("toggle-dark-mode", "value"),
                Input("tertiary-tabs", "active_tab"),
                Input("gsva-stat_df_grid", "columnState"),
                Input("gsva-stat_df_grid", "filterModel"),
                prevent_initial_call=True
            )
            def gsva_stat_df_grid_callback(
                df,
                pathway,
                light,
                active_tab,
                column_state,
                filter_model
            ):

                # import data
                df = pl.LazyFrame.deserialize(BytesIO(b64decode(df)))

                # clean data
                df = df.sort(
                    "Pathway Label"
                ).filter(
                    pl.col("Pathway Label").is_not_null()
                )

                # table theme
                className = "ag-theme-alpine" if light else "ag-theme-alpine-dark"

                if pathway is None:
                    scroll_index = None
                    row_style = None
                else:
                    # filter for feature
                    filtered_df = df.with_row_index().filter(
                        pl.col("Pathway") == pathway
                    ).collect()

                    if filtered_df.is_empty():
                        scroll_index = None
                        row_style = None
                    else:

                        # make scroll index
                        data_row = df.filter(
                            pl.col("Pathway") == pathway
                        )

                        scroll_index = {
                            'data': data_row.collect().to_dicts()[0],
                            'rowPosition': 'middle'
                        }

                        # make row style
                        row_style = {
                            "styleConditions": [{
                                "condition": f"params.data.Pathway === '{pathway}'",
                                "style": {"backgroundColor": "rgba(25, 118, 210, 0.40)"}
                            }]
                        }

                # sort and define all column fields
                df = df.select(
                    ["Pathway", "Pathway Label", pl.all().exclude(["Pathway", "Pathway Label"])]
                )

                fields = make_ag_grid_column_filters(df, exclude=["Pathway"])

                # trigger callback if active tab
                if active_tab == "Pathway Table":
                    return df.collect().to_dicts(), fields, className, row_style, scroll_index

                # return tbl properties
                return df.collect().to_dicts(), fields, className, row_style, scroll_index

        pathway_selection_mediator = adm.PathwaySelectionMediator()
        pathway_selection_mediator.mediate(
            gsea_plot_id=volcano_gsea.id,
            gsea_grid_id='gsea_df_grid',
            gsva_plot_id=volcano_gsva.id,
            gsva_grid_id='gsva-stat_df_grid',
            dropdown_covariates=components["dropdown_covariates"],
            dropdown_pathway=components["dropdown_pathway"]
        )


        if components.get("boxplot_gsva"):

            boxplot_gsva = components["boxplot_gsva"]

            @callback(
                Output(boxplot_gsva.id + '_container', "children"),
                Input("gsva-score-df-store", "data"),
                Input("pathway-selection-store", "data"),
                State('pathway-label-store', 'data'),
                Input("toggle-dark-mode", "value"),
                prevent_initial_call=True
            )
            def plot_pb(df, pathway, pathway_label_store, light):

                # import data
                df = deserialize_lazyframe(df)

                # cast float64 to float32 then actualize DF
                df = df.cast({pl.Float64: pl.Float32}).collect()

                # Hide the plot if a pathway isn't selected.
                if pathway is None:
                    out = [
                        html.Div('Select a Pathway to view GSVA Scores.'),
                        html.Div(style = {'display': 'none'}, children = dcc.Graph(id = boxplot_gsva.id))
                    ]

                    return out

                # de-serialize pathway label store
                pathway_label_store = loads(b64decode(pathway_label_store.encode('utf-8')).decode('utf-8'))
                pathway_label = pathway_label_store[pathway]

                # create title
                #pathway=df['Pathway'].unique()[0]
                plot_title=f"Selected Pathway: <b>{pathway_label}</b>"

                # define templates
                template="plotly_white" if light else "plotly_dark"

                # Group by donor and collect their cytokine values
                # TODO: optimize plotting strategy
                fig = go.Figure()

                fig = boxplot_gsva.plot(
                    df=df,
                    title=plot_title,
                    template=template,
                    condition="Visit",
                    perturbation=False,
                    perturbation_subject_label="Subject ID",
                    perturbation_control_label=config.experiment.control
                )

                out = [dcc.Graph(
                    id = boxplot_gsva.id,
                    figure = fig,
                    style = {
                        'height': "100%",
                        'width': "100%"
                    },
                    config={
                        'displaylogo': False
                    }
                )]

                return out

            @callback(
                Output("gsva-score_df_grid", "rowData"),
                Output("gsva-score_df_grid", "columnDefs"),
                Output("gsva-score_df_grid", "className"),
                Output("gsva-score_df_grid", "getRowStyle"),
                Output("gsva-score_df_grid", "scrollTo"),
                Input("gsva-score-df-store", "data"),
                Input("barcode", "data"),
                Input("toggle-dark-mode", "value"),
                Input("tertiary-tabs", "active_tab"),
                prevent_initial_call=True
            )
            def pb_df_grid_callback(
                df,
                barcode,
                light,
                active_tab
            ):

                # import data
                df = pl.LazyFrame.deserialize(BytesIO(b64decode(df)))

                # define all column fields
                fields = make_ag_grid_column_filters(df, exclude="barcodes")

                # no update if barcode is not formed well
                if barcode is None:
                    return no_update

                # table theme
                className = "ag-theme-alpine" if light else "ag-theme-alpine-dark"

                # make row style
                row_style = None

                # adjust tbl scroll
                scroll_index = None

                # trigger callback if active tab
                if active_tab == "Pseudo-Bulk Table":
                    return df.collect().to_dicts(), fields, className, row_style, scroll_index

                # return tbl properties
                return df.collect().to_dicts(), fields, className, row_style, scroll_index


        if config.data.deg.specifications.deg_by_covariates:

            deg_by_covariates_plot = components["deg_by_covariates_plot"]
            deg_by_covariates_modal = components["deg_by_covariates_modal"]
            deg_by_covariates_button = components["deg_by_covariates_button"]

            # turn on modal given a click of corresponding button
            @callback(
                Output(deg_by_covariates_modal.id, 'is_open'),
                Input(deg_by_covariates_button.id, 'n_clicks'),
                prevent_initial_call=True
            )
            def show_gex_by_all_conditions(button):
                if button:
                    return True
                return False

            # make the plot
            @callback(
                Output(deg_by_covariates_plot.id, "figure"),
                Input("deg-by-covariates-df-store", "data"),
                Input("toggle-dark-mode", "value"),
                Input("gene-selection-store", "data"),
                Input("covariate-labels", "data"),
                Input("covariate-values", "data"),
                prevent_initial_call=True
            )
            def plot_deg_by_covariates_plot(df, light, gene, labels, values):

                # import data
                df = deserialize_lazyframe(df).collect()

                template = "plotly_white" if light else "plotly_dark"

                if gene is None:
                    gene="PDE4B"

                # return plot with error message if gene wasn't found.
                if df[deg_by_covariates_plot.nlp_col].is_null().all():
                    if df[deg_by_covariates_plot.es_col].is_null().all():
                        fig = scatter(data_frame=pl.DataFrame({}), template=template)
                        fig.add_annotation(
                            x=0.50,
                            y=0.25,
                            xref='x domain',
                            yref='y domain',
                            text=f"Gene expression results for {gene} could not <br>be found across experimental conditions.",
                            font=dict(
                                size=18,
                                weight="bold"
                            ),
                            showarrow=False
                        )
                        return fig

                fig = deg_by_covariates_plot.plot(
                    data_frame=df,
                    template=template,
                    custom_data=["Label"]
                )

                fig.update_traces(
                    hovertemplate=
                        '<b>%{customdata[0]}</b>' + '<br><br>' +
                        f"{deg_by_covariates_plot.nlp_col} = " + "%{y:.2f}<br>" +
                        f"{deg_by_covariates_plot.es_col} = " + "%{x:.3f}<br>" +
                        '<extra></extra>'
                )

                # annotate the current/selected covariate values
                if values:
                    values = [ i for i in values if "censor" not in i.lower() ]
                    highlight_point_with_label(
                        df = df,
                        fig = fig,
                        feature_col="Label",
                        feature_value=": ".join(values),
                        x_col=deg_by_covariates_plot.es_col,
                        y_col=deg_by_covariates_plot.nlp_col,
                        color_map=deg_by_covariates_plot.point_colors,
                        title=True
                    )
                
                # add gene expression direction labels
                fig.add_annotation(text=f'{gene} higher in Foreground', showarrow=False, x=1, y=1.05, xanchor='right', yanchor='top', xref='paper', yref='paper', font=dict(size=16))
                fig.add_annotation(text=f'{gene} higher in Background', showarrow=False, x=0, y=1.05, xanchor='left', yanchor='top', xref='paper', yref='paper', font=dict(size=16))

                return fig

            # allow click to change dropdown covariates. Acts like a mediator!
            @callback(
                [Output(comp.id, "value", allow_duplicate=True) for comp in components["dropdown_covariates"]],
                State("covariate-values", "data"),
                Input(deg_by_covariates_plot.id, "clickData"),
                prevent_initial_call=True
            )
            def deg_by_covariates_mediate(values, click):
                if click:
                    covariates = click["points"][0]["customdata"][0].split(": ")
                    # covariates.extend([values[-1]])
                    return covariates

        if config.data.deg.specifications.summaries_plot:

            deg_summaries_plot = components["deg_summaries_plot"]
            deg_summaries_modal = components["deg_summaries_modal"]
            deg_summaries_button = components["deg_summaries_button"]
            deg_summaries_store = components["deg-summaries-store"]

            @callback(
                Output(deg_summaries_modal.id, 'is_open'),
                Input(deg_summaries_button.id, 'n_clicks'),
                prevent_initial_call=True
            )
            def show_deg_summaries_modal(button):
                if button:
                    return True
                return False

            @callback(
                Output("deg-summaries-store", "data"),
                Input("covariate-values", "data"),
                Input("covariate-labels", "data"),
                Input("barcode", "data"),
                prevent_initial_call=True
            )
            def generate_deg_summaries_df(values, labels, barcode):

                if barcode is None:
                    return no_update

                # load data
                df = pl.scan_parquet(config.data.deg.files.summaries)

                # get DF using covariates that don't vary like cytokine or treatment.
                for covariate in deg_summaries_plot.const_covariate:

                    # find corresponding value of the covariate we want to make constant
                    v = values[labels.index(covariate)]

                    df = df.filter(pl.col(covariate) == v)

                # wide to long format
                df = df.unpivot(index=["Cytokine", "Cell Type", "Censoring", "N Tested", "Type Color", "N Up", "N Down"])

                df = df.cast({pl.UInt32: pl.Int32})

                # return data
                serialized_df = serialize_lazyframe(df)

                return serialized_df

            @callback(
                Output(deg_summaries_plot.id, "figure"),
                Input(deg_summaries_store.id, "data"),
                Input("covariate-values", "data"),
                Input("covariate-labels", "data"),
                Input("toggle-dark-mode", "value"),
                prevent_initial_call=True
            )
            def make_deg_summaries_plot(df, values, labels, light):

                # import data
                df = deserialize_lazyframe(df).collect()

                # define templates
                template="plotly_white" if light else "plotly_dark"

                # split df into up and down regulated genes
                up_df = df.select(
                    labels + ["N Up", "Type Color"]
                )

                dn_df = df.select(
                    labels + ["N Down", "Type Color"]
                )

                # set down regulated gene colors to half opacity
                dn_df = dn_df.with_columns(
                    pl.col("Type Color").map_elements(
                        hex_to_rgba, return_dtype=pl.String
                    ).alias("Type Color")
                )

                # plot up-regulated genes first
                fig = go.Figure(
                    data= [
                        go.Bar(
                            x=up_df[deg_summaries_plot.vary_covariates[0]],
                            y=up_df["N Up"],
                            marker_color=df["Type Color"],
                            text=up_df["N Up"],
                            textposition="outside",
                            hovertemplate='%{x}' + '<extra></extra>'
                        )
                    ]
                )

                # then add down-regulated genes
                fig.add_trace(
                    go.Bar(
                        x=dn_df[deg_summaries_plot.vary_covariates[0]],
                        y=-1*dn_df["N Down"],
                        marker_color=dn_df["Type Color"],
                        text=dn_df["N Down"],
                        textposition="outside",
                        hovertemplate='%{x}' + '<extra></extra>'
                    )
                )

                fig.update_traces(
                    hoverlabel = {
                        'align': 'left',
                        'bordercolor': 'rgb(0,0,0,0)',
                        'font': {
                            'color': 'white',
                            'size': 16,
                            'weight': 600,
                            'family': 'Roboto Condensed'
                        }
                    }
                )

                # add up and down-regulated genes label
                fig.add_annotation(
                    text=f"Higher in {values[0]}",
                    xref="paper",
                    yref="paper",
                    x=0.5,
                    y=1,
                    xanchor="center",
                    yanchor="top",
                    showarrow=False,
                    font=dict(
                        size=16
                    )
                )

                fig.add_annotation(
                    text="Higher in PBS",
                    xref="paper",
                    yref="paper",
                    x=0.5,
                    y=0,
                    xanchor="center",
                    yanchor="bottom",
                    showarrow=False,
                    font=dict(
                        size=16
                    )
                )

                fig.update_xaxes(
                    tickangle=-45,
                    automargin=True
                )

                fig.update_yaxes(
                    automargin=True
                )

                fig.update_layout(
                    barmode="relative",
                    hoverlabel=dict(bordercolor="rgb(0,0,0,0)"),
                    clickmode = 'event',
                    autosize=True,

                    margin=dict(t=0,b=0,l=0,r=0),

                    font_family='Roboto Condensed',
                    font_size=14,

                    template = template,
                    showlegend=False
                )

                return fig

            # the following callback also acts like a mediator.
            @callback(
                [Output(comp.id, "value", allow_duplicate=True) for comp in components["dropdown_covariates"]],
                State("covariate-values", "data"),
                Input(deg_summaries_plot.id, "clickData"),
                prevent_initial_call=True
            )
            def DegSummariesCovariatesMediator(values, click):
                if click:
                    values[1] = click['points'][0]['label']
                    return values

    # create callbacks for the layout
    create_layout_callbacks()

    # create callbacks for minimal deg functionality
    create_minimal_deg_callbacks(components=components)

    # create callbacks for optional deg functionality
    create_optional_deg_callbacks(config, components=components)

    return None
