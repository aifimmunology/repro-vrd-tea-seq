import polars as pl
from plotly.express import scatter
from plotly.graph_objects import Figure
from dash_bootstrap_components import Card, CardHeader, CardBody
from dash.dcc import Graph, Loading

class MAPlot:
    def __init__(
        self,
        id:             str,
        section_title:  str,
        means_col:      str,
        es_col:         str
    ):
        self.id = id
        self.section_title = section_title
        self.means_col = means_col
        self.es_col = es_col

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
        es_cutoff:      float,
        custom_data:    list,
        **kwargs
    ) -> Figure:

        fig = scatter(
            data_frame=data_frame,
            x=self.means_col,
            y=self.es_col,
            color='Category',
            color_discrete_map = self.point_colors,
            render_mode='webgl',
            custom_data=custom_data
        )

        # add horizontal enrichment score cutoff lines
        fig.add_hline(y=es_cutoff, line_dash="dash", line_color='rgba(148, 148, 148, 0.58)')
        fig.add_hline(y=-1*es_cutoff, line_dash="dash", line_color='rgba(148, 148, 148, 0.58)')

        fig.update_traces(
            hovertemplate=
                '<b>%{customdata[0]}</b>' + '<br><br>' +
                '<b>Mean Expression</b>: %{customdata[1]:.2f}' + '<br>' +
                '<b>Log2FC</b>: %{customdata[2]:.2f}' + '<br>' +
                "<extra></extra>",
            hoverlabel = {
                'align': 'left',
                'font': {
                    'color': 'white',
                    'size': 14,
                    'weight': 600,
                }
            }
        )

        fig.update_xaxes(
            range=[
                data_frame.select(pl.min(self.means_col))[0] * 0.90,
                data_frame.select(pl.max(self.means_col))[0] * 1.10,
            ],
            automargin=True
        )

        fig.update_yaxes(
            automargin=True
        )

        fig.update_layout(
            hoverlabel=dict(bordercolor="rgb(0,0,0,0)"),
            clickmode = 'event',
            font_family='Roboto Condensed',
            font_size=14,

            autosize = True,

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

    def build_ui(self):
        ui = Card([
            CardHeader(self.section_title),
            CardBody([
                Loading(
                    Graph(
                        id = self.id,
                        figure = self.fig,
                        style = {
                            'height': '100%',
                            'width': '100%'
                        }
                    ),
                overlay_style={
                    "visibility":"visible",
                    "filter": "blur(2px)"
                },
                type='circle',
                delay_show=2000
                ) # end Loading
            ]) # end CardBody
        ],
        style={
            'height': '100%',
            'width': '100%'
        } # end Card children
    ) # end Card

        self.ui = ui
