import polars as pl
from plotly.express import scatter
from plotly.graph_objects import Figure
from dash_bootstrap_components import Card, CardHeader, CardBody
from dash.dcc import Graph, Loading

class VolcanoPlot:
    def __init__(
        self,
        id:             str,
        section_title:  str,
        p_col:          str,
        es_col:         str
    ):
        self.id = id
        self.section_title = section_title
        self.p_col = p_col
        self.es_col = es_col
        self.nlp_col = f"-log10({self.p_col})"

        self.fig = {}

        self.point_colors:  dict = {
            'down-regulated': '#3577B4',
            'up-regulated': '#EB6648',
            'not significant': '#828282'
        }

        self.build_ui()

    def plot(
        self,
        data_frame:     pl.DataFrame,
        template:       str,
        custom_data:    list,
        **kwargs
    ) -> Figure:

        # Generate plot
        fig = scatter(
            data_frame = data_frame,
            x = self.es_col,
            y = self.nlp_col,
            color = 'Category',
            color_discrete_map = self.point_colors,
            render_mode = 'webgl',
            custom_data=custom_data
        )

        # add vertical and horizontal lines along axes for styling
        fig.add_vline(x=0, line_width=1, opacity=0.50)
        fig.add_hline(y=0, line_width=1, opacity=0.50)

        # set x and y axes ranges
        fig.update_xaxes(
            range=self.determine_axis_range(data_frame, self.es_col, type='xaxis'),
            automargin=True
        )
        fig.update_yaxes(
            range=self.determine_axis_range(data_frame, self.nlp_col, type='yaxis'),
            automargin=True
        )

        # add hover template and labels if provided
        try:
            if kwargs["hovertemplate"]:
                fig.update_traces(
                    hovertemplate=kwargs["hovertemplate"]
                )
        except KeyError:
            pass

        fig.update_traces(
            hoverlabel = {
                'align': 'left',
                'bordercolor': 'rgb(0,0,0,0)',
                'font': {
                    'color': 'white',
                    'size': 14,
                    'weight': 600,
                    'family': 'Roboto Condensed'
                }
            }
        )

        # additional styling
        fig.update_layout(
            clickmode = 'event',
            font_family='Roboto Condensed',
            font_size=14,

            autosize = True,
            # margin=dict(t=75,l=75,r=75,b=75),

            xaxis=dict(title_standoff=10),
            yaxis=dict(title_standoff=25),

            template = template,
            legend = dict(
                itemsizing='constant',
                orientation='h',
                xref='container',
                yref='container',
                yanchor="bottom",
                xanchor="center",
                y=-0.25,
                x=0.50
            )
        )

        return fig

    def determine_axis_range(self, data: pl.DataFrame, column: str, type="xaxis"):

        max_value = data.select(pl.col(column).abs().max())[0, 0]

        if max_value is not None:
            max_value = max_value * 1.10
        else:
            max_value = 2

        if type == 'xaxis':
            return [-max_value, max_value]
        else:
            return [0, max_value]

    def build_ui(self):
        ui = Card([
            CardHeader(self.section_title),
            CardBody([
                Loading(
                    Graph(
                        id = self.id,
                        figure = self.fig,
                        style = {
                            'height': "100%",
                            'width': "100%"
                        },
                        config={
                            'displaylogo': False
                        }
                    ),
                overlay_style={
                    "visibility":"visible",
                    "filter": "blur(2px)"
                },
                type='circle',
                delay_show=1500
                ) # end Loading
            ]) # end CardBody
        ],
        style={
            'height': '100%',
            'width': '100%'
        } # end Card children
    ) # end Card

        self.ui = ui
