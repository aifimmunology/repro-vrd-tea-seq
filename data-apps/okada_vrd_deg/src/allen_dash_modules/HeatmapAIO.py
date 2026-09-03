import uuid

import numpy as np
import polars as pl
import plotly.express as px

import dash_bootstrap_components as dbc

from dash_iconify import DashIconify

from dash import Input, Output, State, callback, html, dcc, MATCH, ctx, Patch

TAB_CONTENT_STYLE = {'padding': '1rem'}
ICON_STYLE = {
    'backgroundColor': '#f8f9fa',
    'margin-right': '0.5rem',
    'margin-left': '0.5rem'
}
CONTROL_COMPONENT_STYLE = {
    'border': 'var(--bs-border-color-translucent)',
    'border-style': 'solid',
    'border-width': '0.1rem',
    'border-radius': 'var(--bs-border-radius)',
    'backgroundColor': '#f8f9fa',
    'padding': '0.5rem',
    'width': '100%',
    'height': '100%'
}

PALETTE_ICON = DashIconify(icon="material-symbols:palette-outline", width=32, height=32, style=ICON_STYLE)
EYEDROPPER_ICON = DashIconify(icon="ph:eyedropper", width=32, height=32, style=ICON_STYLE)

class HeatmapAIO(html.Div):

    class ids:

        df = lambda aio_id: {
            'component': 'HeatmapAIO',
            'subcomponent': 'df',
            'aio_id': aio_id
        }

        tabs = lambda aio_id: {
            'component': 'HeatmapAIO',
            'subcomponent': 'tabs',
            'aio_id': aio_id
        }

        graph_display = lambda aio_id: {
            'component': 'HeatmapAIO',
            'subcomponent': 'graph-display',
            'aio_id': aio_id
        }

        store = lambda aio_id: {
            'component': 'HeatmapAIO',
            'subcomponent': 'store',
            'aio_id': aio_id
        }

        example_display = lambda aio_id: {
            'component': 'HeatmapAIO',
            'subcomponent': 'example-display',
            'aio_id': aio_id
        }

        example_button = lambda aio_id: {
            'component': 'HeatmapAIO',
            'subcomponent': 'example-button',
            'aio_id': aio_id
        }

        color_dropdown = lambda aio_id: {
            'component': 'HeatmapAIO',
            'subcomponent': 'color-dropdown',
            'aio_id': aio_id
        }

        palette_type = lambda aio_id: {
            'component': 'HeatmapAIO',
            'subcomponent': 'palette-type',
            'aio_id': aio_id
        }

    ids = ids

    def __init__(
        self,
        df: pl.DataFrame | pl.LazyFrame | str = None,
        aio_id: str = None,
        annotation_cols: list = None,
        **kwargs
    ):

        # 1. initiate aio_id
        if not aio_id:
            aio_id = str(uuid.uuid4())

        # 2. Define data store
        if isinstance(df, pl.LazyFrame):
            dataframe = df.collect().to_dicts()
        elif isinstance(df, pl.DataFrame):
            dataframe = df.to_dicts()
        else:
            raise Exception("Input dataframe is neither polars LazyFrame or DataFrame!")

        store_component = dcc.Store(data=dataframe, id=self.ids.store(aio_id=aio_id))

        # 3. Define layout
        data_tab = dbc.Tab(
            children=[
                html.P("Content of tab 1"),
                dbc.Button("Randomize Data!", id=self.ids.example_button(aio_id=aio_id))
            ],
            style=TAB_CONTENT_STYLE
        )

        palette_type_selector = html.Div([
            html.I(PALETTE_ICON),
            html.Div(dbc.RadioItems(
                id=self.ids.palette_type(aio_id=aio_id),
                options=["Diverging", "Sequential"],
                value="Diverging",
                inline=True,
                style={'marginTop': '0.25rem', 'marginLeft': '0.5rem'}
            ), style={'width': '100%'})
        ],
            style={'display': 'inline-flex'} | CONTROL_COMPONENT_STYLE,
            className='dbc'
        )

        # Get initial diverging options
        initial_options = sorted([
            k for k in dir(px.colors.diverging)
            if not k.startswith("_") and isinstance(getattr(px.colors.diverging, k), list)
        ])

        color_palette_selector = html.Div([
            html.I(EYEDROPPER_ICON),
            html.Div(dcc.Dropdown(
                id=self.ids.color_dropdown(aio_id=aio_id),
                options=initial_options,
                multi=False,
                clearable=False,
                value=initial_options[0]
            ), style={'width': '100%'})
        ],
            style={'display': 'inline-flex'} | CONTROL_COMPONENT_STYLE,
            className='dbc'
        )

        color_tab = dbc.Tab([
            dbc.Row([
                dbc.Col(palette_type_selector, width=6, xs=12, sm=12, md=12, lg=12, xl=6, xxl=4),
                dbc.Col(color_palette_selector, width=6, xs=12, sm=12, md=12, lg=12, xl=6, xxl=8)
            ])
        ], style=TAB_CONTENT_STYLE)

        settings_tab = dbc.Tab(children=["Content of tab 3"], style=TAB_CONTENT_STYLE)

        tab_set = dbc.Tabs(
            id=self.ids.tabs(aio_id=aio_id),
            children=[
                data_tab,
                color_tab,
                settings_tab
            ],
            className='heatmap-aio-tabs nav-pills nav-fill'
        )

        left_pane = html.Div([
            tab_set
        ], style={
            'width': '40%',
            'margin': '0.5rem',
            'borderWidth': '0.1rem',
            'borderRadius': '0.5rem',
            'borderColor': '#808080',
            'borderStyle': 'dashed'
        })

        right_pane = html.Div([
            dcc.Loading(
                dcc.Graph(
                    id=self.ids.graph_display(aio_id=aio_id),
                    figure=None
                ),
            delay_show=1500,
            type="circle")
        ], style={
            'width': '60%',
            'margin': '0.5rem',
            'borderWidth': '0.1rem',
            'borderRadius': '0.5rem',
            'borderColor': '#808080',
            'borderStyle': 'dashed'
        })

        bottom_pane = html.Div([
            "Default value"
        ], id=self.ids.example_display(aio_id=aio_id))

        # 4. Construct component layout
        ui = html.Div([
            store_component,
            html.B(f"Heatmap {aio_id}"),
            html.Div([left_pane, right_pane], style={'display': 'inline-flex', 'width': '100%'}),
            bottom_pane
        ],
        style={
            'margin': '1rem',
            'padding': '1rem',
            'borderWidth': '0.1rem',
            'borderRadius': '0.5rem',
            'borderColor': '#808080',
            'borderStyle': 'solid'
            }
        )

        self.ui = ui

    # callbacks ---------------------------------------------------------------

    # example callback for example_button
    @callback(
        Output(ids.example_display(MATCH), 'children'),
        Input(ids.example_button(MATCH), 'n_clicks'),
        prevent_initial_call=True
    )
    def example_callback(clicks):
        return f"You clicked this button {clicks} times!"

    # change color palette options given type
    @callback(
        Output(ids.color_dropdown(MATCH), 'options'),
        Output(ids.color_dropdown(MATCH), 'value'),
        Input(ids.palette_type(MATCH), 'value'),
        prevent_initial_call=True
    )
    def switch_color_dropdown_choices(value):
        if value == "Sequential":
            source = px.colors.sequential
        else:
            source = px.colors.diverging

        options = sorted([k for k in dir(source) if not k.startswith("_") and isinstance(getattr(source, k), list)])
        return options, options[0]

    # plotting the heatmap
    @callback(
        Output(ids.graph_display(MATCH), 'figure'),
        State(ids.store(MATCH), 'data'),
        Input(ids.color_dropdown(MATCH), 'value'),
        Input(ids.example_button(MATCH), 'n_clicks'),
        State(ids.palette_type(MATCH), 'value'),
        prevent_initial_call="initial_duplicate"
    )
    def plot_heatmap(df, palette, clicks, palette_type):

        # import data
        df = pl.DataFrame(df)

        # random sample rows
        sample_indices = np.random.choice(range(0, df.shape[0]), size=10, replace=False)
        df = df.with_row_index().filter(pl.col("index").is_in(sample_indices)).drop("index")

        # Resolve color scale
        if palette_type == "Sequential":
            color_scale = getattr(px.colors.sequential, palette)
        else:
            color_scale = getattr(px.colors.diverging, palette)

        # construct figure
        fig = px.imshow(df, color_continuous_scale=color_scale)

        return fig
