from base64 import b64encode

from dash import html, dcc

import dash_bootstrap_components as dbc

from dash_iconify import DashIconify

import polars as pl

import allen_dash_modules as adm

def create_deg_layout(config: adm.ConfigSingleton, components: dict) -> list:

    def make_stores(config: adm.ConfigSingleton) -> list:
        """
        These are the minimum set of stores to make the DEG work. Additional
        stores will be added optionally as specified by an input JSON config file.
        """

        # read and encode DEG meta.parquet file
        deg_meta_parquet = pl.scan_parquet(config.data.deg.files.meta)
        deg_meta_parquet = b64encode(deg_meta_parquet.serialize()).decode('utf-8')

        return [
            dcc.Store("volcano-deg-df-1"),
            dcc.Store("volcano-deg-df"),

            dcc.Store("p-cutoff-store"),
            dcc.Store("es-cutoff-store"),
            dcc.Store("gene-selection-store", data=None),

            dcc.Store("covariates", data=None),
            dcc.Store("covariate-values"),
            dcc.Store("covariate-labels", data=config.experiment.covariates),
            dcc.Store("barcode"),
            dcc.Store("covariates-df-store", data=deg_meta_parquet),

            dcc.Store("pathway-selection-store", data='INFLAMMATORY RESPONSE'),
            dcc.Store("pathway-label-store", data=None),
            dcc.Store("pathway-gene-store", data=None),

            dcc.Store("startup-modal-store", storage_type='local', data=True)
        ]

    def make_header(config: adm.ConfigSingleton) -> list:

        """
        Make a header using pre-defined ADM components.
        Consider using flags to set theme, setting, and help in the input JSON + fx args.
        """

        # basic visual components
        header = adm.AppHeader(
            toggle_theme=True,
            toggle_settings=True,
            toggle_help=True
        )

        header.ui = header.build_app_header(
            title=config.dashboard.project_title,
            institute_url=config.dashboard.institute_homepage,
            study_url=config.dashboard.project_homepage
        )

        # instantiate help page with a markdown file
        header.help_page.ui = header.help_page.build_help_panel(
            markdown_file=config.dashboard.help_page_content
        )

        # instantiate components to append to settings page
        # for example, add a checkbox to enable the sticky pathway feature
        settings_page_content = [html.H2("Plot Options")]
        # if 'gsea_rank_toggle' in components.keys():
        #     settings_page_content.extend(components["gsea_rank_toggle"])
        settings_page_content.extend(components["sticky_pathway_checkbox"])


        header.settings_page.ui = header.settings_page.build_settings_panel(
            children=dbc.Offcanvas(
                children=settings_page_content
            )
        )

        return [header.ui, header.settings_page.ui, header.help_page.ui]

    def make_layout(config: adm.ConfigSingleton, components: dict) -> list:

        """
        Create minimal components like dropdowns, plots, and tables for the app to work.
        """

        def make_startup_modal(config: adm.ConfigSingleton) -> adm.StartupModal:

            if config.dashboard.startup_modal:
                startup_modal = adm.StartupModal()
                startup_modal.ui = startup_modal.build_startup_modal(
                    markdown_file=config.dashboard.startup_modal_content
                )
                return startup_modal
            else:
                return []

        def make_sidebar(components: dict) -> html.Div:
            """
            Create the sidebar with optional components
            """

            close_sidebar_icon = DashIconify(icon="bxs:left-arrow")

            sidebar_headers = {
                "header1": html.Div([
                    html.H5(
                        [
                            "Experimental Conditions",
                            DashIconify(
                                icon = "material-symbols:help-outline",
                                id="experimental-conditions-popover",
                                className='iconify-in-text'
                            )
                        ]
                    ),
                    dbc.Popover(
                        [
                            dbc.PopoverBody(
                                [
                                    html.Span("This section adjusts the major experimental conditions that influence several downstream plots.\n\n"),
                                    html.Span([
                                        html.B("Comparison type: "), html.Span("make a cross-sectional (between cohorts within the same timepoint) or longitudinal (same cohort between timepoints) comparison.\n\n")
                                    ]),
                                    html.Span([
                                        html.B("Foreground: "), html.Span("select the condition at the right side of the volcano plot.\n\n")
                                    ]),
                                    html.Span([
                                        html.B("Background: "), html.Span("select the condition at the left side of the volcano plot.\n\n")
                                    ]),
                                    html.Span([
                                        html.B("Cell Type: "), html.Span("subset results by cell type.\n\n")
                                    ])
                                ]
                            )
                        ],
                        target="experimental-conditions-popover",
                        trigger="hover"
                    )
                ]),
                
                "header2": html.Div([
                    html.H5(
                        [
                            "Significance Cutoffs",
                            DashIconify(
                                icon = "material-symbols:help-outline",
                                id="significance-cutoffs-popover",
                                className='iconify-in-text'
                            )
                        ]
                    ),
                    dbc.Popover(
                        [
                            dbc.PopoverBody("Adjust the p-value and the log2fc values to label genes as up- or down-regulated."
                            )
                        ],
                        target="significance-cutoffs-popover",
                        trigger="hover"
                    )
                ]),

                "header3": html.Div([
                    html.H5(
                        [
                            "Pathways",
                            DashIconify(
                                icon = "material-symbols:help-outline",
                                id="pathway-popover",
                                className='iconify-in-text'
                            )
                        ]
                    ),
                    dbc.Popover(
                        [
                            dbc.PopoverBody("Select the source of a GSEA gene set (either Hallmark or \
                                            Reactome) and a pathway using the following dropdown OR under \
                                            the 'Pathway Label' column under the 'Pathway Table' tab.")
                        ],
                        target="pathway-popover",
                        trigger="hover"
                    )
                ]),

                "header4": html.Div([
                    html.H5(
                        [
                            "Gene",
                            DashIconify(
                                icon = "material-symbols:help-outline",
                                id="genes-popover",
                                className='iconify-in-text'
                            )
                        ],
                    ),
                    dbc.Popover(
                        [
                            dbc.PopoverBody("Show a gene's expression via the following dropdown OR the 'Gene' column under the 'DEG Table' tab.")
                        ],
                        target="genes-popover",
                        trigger="hover"
                    )
                ])

            }

            # assemble minimal sub-sections of the sidebar
            covariates_section = [sidebar_headers["header1"]] + [ comp.ui for comp in components["dropdown_covariates"]] + [html.Hr(style={'margin': '1rem 0rem 1rem 0rem'})]
            cutoff_section = [sidebar_headers["header2"]] + [ comp.ui for comp in components["dropdown_plot_filters"] ] + [html.Hr(style={'margin': '1rem 0rem 1rem 0rem'})]
            gene_section = [sidebar_headers["header4"], components["dropdown_gene"].ui]

            # insert additional logic here to supplement sidebar sub-sections with optional components
            if components.get("deg_by_covariates_plot"):
                gene_section.append(components["deg_by_covariates_button"])

            if components.get("deg_summaries_plot"):
                gene_section.append(components["deg_summaries_button"])
            
            gene_section.append(html.Hr())

            if components.get("dropdown_pathway_source"):
                pathway_section = [
                    sidebar_headers["header3"],
                    components["dropdown_pathway_source"].ui,
                    components["dropdown_pathway"].ui,
                ] + components["gsea_rank_toggle"]
            else:
                pathway_section = []

            # create the sidebar
            sidebar = html.Div([

                html.Div(children=[

                    # section 1: experimental conditions
                    html.Div(covariates_section, className='sidebar-section'),

                    # section 2: pval and log2fc score cutoffs
                    html.Div(children=cutoff_section, className='sidebar-section'),

                    # section 4: gene
                    html.Div(children=gene_section, className='sidebar-section'),

                    # section 3: pathway
                    html.Div(children=pathway_section, className='sidebar-section'),

                ], id="sidebar", className="sidebar"),

                # pull tab OUTSIDE the sidebar
                html.Div(close_sidebar_icon, id="pull-tab", className="pull-tab"),
                dbc.Popover(
                    "Collapse side panel",
                    target="pull-tab",
                    id="pull-tab-popover",
                    trigger="hover",
                    body=True
                )

            ])

            return sidebar

        def make_windows(components: dict) -> html.Div:

            """
            This function creates all the individual windows and then
            assembles all the plots together with the correct styling.
            """

            # define primary window given components
            primary_window = components["volcano_deg"].ui
            secondary_window = []
            secondary_window_content = {}
            tertiary_window = []
            tertiary_window_content = {}

            # define secondary windows given components...
            if components.get("volcano_gsea"):
                secondary_window_content["GSEA Stats"] = components["volcano_gsea"].ui
            
            if components.get("volcano_gsva"):
                secondary_window_content["GSVA Stats"] = components["volcano_gsva"].ui
            
            if components.get("boxplot_gsva"):
                secondary_window_content["GSVA Scores"] = components["boxplot_gsva"].ui
            
            if components.get("boxplot_pb"):
                secondary_window_content["Pseudo-Bulk"] = components["boxplot_pb"].ui
            
            if components.get("ma_plot"):
                secondary_window_content["MA Plot"] = components["ma_plot"].ui
            
            secondary_window = create_tabbed_windows(
                secondary_window_content,
                tabset_id='secondary-tabs'
            )

            # define tertiary windows given components...
            if components.get("gene_panel"):
                tertiary_window_content["Gene Info"] = html.Div(id=components["gene_panel"].output_div)

            tertiary_window_content["DEG Table"] = html.Div(components["volcano_deg_grid"])

            if components.get("volcano_gsea"):
                tertiary_window_content["GSEA Results Table"] = html.Div(components["gsea_df_grid"])
            
            if components.get("volcano_gsva"):
                tertiary_window_content["GSVA Results Table"] = html.Div(components["gsva-stat_df_grid"])
            
            if components.get("boxplot_gsva"):
                tertiary_window_content["GSVA Scores Table"] = html.Div(components["gsva-score_df_grid"])

            if components.get("boxplot_pb"):
                tertiary_window_content["Pseudo-Bulk Table"] = html.Div(components["pb_df_grid"])

            tertiary_window_content["Covariates"] = html.Div(components["covariates_df_grid"])

            # place tertiary window contents into a tab layout.
            tertiary_window = create_tabbed_windows(
                components=tertiary_window_content,
                tabset_id='tertiary-tabs'
            )

            # assemble primary and secondary window into top_windows
            top_windows = html.Div(children=[
                html.Div(primary_window, className='primary-window'),
                html.Div(secondary_window, className='secondary-window')
            ], className='top-windows')

            # define bottom windows aka html.Div(tertiary windows) 
            bottom_windows = html.Div(tertiary_window, className='tertiary-window')

            # assemble all windows together into one Div
            all_windows = html.Div([top_windows, bottom_windows], id='main-content', className='main-content')

            return all_windows

        # invoke the following sub-functions to create the app layout!

        # 1. create optional startup modal ------------------------------------
        startup_modal = make_startup_modal(config)

        # 2. create layout components -----------------------------------------
        sidebar = make_sidebar(components=components)
        container = make_windows(components=components)
        barcode_alert = components["barcode_alert"]

        # 3. assemble entire visual layout
        layout = html.Div(children=[sidebar, container], className='container')

        # 4. include all hidden Ag.Grids and optional stores ------------------
        # all_tables = []
        all_optional_stores = []

        for component in components.values():
            if isinstance(component, dcc.Store):
                all_optional_stores.append(component)

        # all_tables = html.Div(all_tables, hidden=True)
        all_optional_stores = html.Div(all_optional_stores, hidden=True)

        # 6. include additional hidden or custom components -------------------
        # these components can be part of the default app view depending on how you make them.
        custom_components = []

        # insert modals (visual bootstrap components, not plots) into custom_components list.
        try:
            if components["deg_by_covariates_modal"]:
                custom_components.append(components["deg_by_covariates_modal"])
        except KeyError:
            pass

        try:
            if components["deg_summaries_modal"]:
                custom_components.append(components["deg_summaries_modal"])
        except KeyError:
            pass

        custom_components = html.Div(custom_components)

        # assemble every component into a list of dbc components or Div -------
        layout = [
            startup_modal.ui,
            barcode_alert,
            layout,
            all_optional_stores,
            custom_components
        ]

        return layout

    # invoke all sub-functions!
    minimal_stores = make_stores(config)
    header = make_header(config)
    layout = make_layout(config, components=components)

    # gather all results into a list of html and dash bootstrap components.
    all_layout_components = minimal_stores + header + layout

    return all_layout_components

# supporting functions ----------------------------------------------------------------------------

def create_tabbed_windows(components: dict, tabset_id: str) -> html.Div:
    """
    Key should be the label of the tab and value of either a Div or dbc component.
    {
        "Label 1": html.Div(id="output-for-data-type-1"),
        "Label 2": components["myplot"].ui,
    }
    """

    if len(components) == 0:
        return []

    tabbed_window = html.Div([
        dbc.Tabs(
            children=[
                dbc.Tab(children=value, label=key) for key, value in components.items()
            ],
            className='nav-pills',
            id=tabset_id
        )
    ])

    return tabbed_window
