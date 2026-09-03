import polars as pl
from plotly.express import bar
from plotly.graph_objects import Figure
from dash_bootstrap_components import Card, CardHeader, CardBody
from dash.html import H4
from dash.dcc import Graph, Loading

class DegSummariesPlot:
    def __init__(
        self,
        id:                 str,
        section_title:      str,
        const_covariates:   list[str],
        vary_covariate:     list[str]
    ):
        self.id = id
        self.section_title = section_title
        self.const_covariates = const_covariates
        self.vary_covariate = vary_covariate

        self.fig = {}

        self.build_ui()

    def plot(
        self,
        data_frame:     pl.DataFrame,
        template:       str,
        **kwargs
    ) -> Figure:

        print(data_frame)

        fig = bar(
            data_frame=data_frame,
            **kwargs
        )

        # fig = scatter(
        #     data_frame=data_frame,
        #     x="means",
        #     y="log2fc",
        #     color='Category',
        #     color_discrete_map = self.point_colors,
        #     render_mode='webgl',
        #     custom_data=custom_data
        # )

        # # add horizontal enrichment score cutoff lines
        # fig.add_hline(y=es_cutoff, line_dash="dash", line_color='rgba(148, 148, 148, 0.58)')
        # fig.add_hline(y=-1*es_cutoff, line_dash="dash", line_color='rgba(148, 148, 148, 0.58)')

        # fig.update_traces(
        #     hovertemplate=
        #         '<b>%{customdata[0]}</b>' + '<br><br>' +
        #         '<b>Mean Expression</b>: %{customdata[1]:.2f}' + '<br>' +
        #         '<b>Log2FC</b>: %{customdata[2]:.2f}' + '<br>' +
        #         "<extra></extra>",
        #     hoverlabel = {
        #         'align': 'left',
        #         'font': {
        #             'color': 'white',
        #             'size': 18,
        #             'weight': 600
        #         }
        #     }
        # )

        # fig.update_xaxes(range=[
        #     data_frame.select(pl.min(self.means_col))[0] * 0.90,
        #     data_frame.select(pl.max(self.means_col))[0] * 1.10,
        # ])

        fig.update_layout(
            hoverlabel = dict(bordercolor="rgb(0,0,0,0)"),
            clickmode = 'event',
            font = {
                    'family': 'Roboto Condensed',
                    'size': 18
            },

            template = template
        )

        return fig

    def build_ui(self):
        ui = Card([
            CardHeader(H4(self.section_title)),
            CardBody([
                Loading(
                    Graph(
                        id = self.id,
                        figure = self.fig,
                        style = {
                            'height': '100%',
                            'width': '100%'
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
                delay_show=2000
                )
            ])
        ],

        style={
            'height': '100%',
            'width': '100%'
        })

        self.ui = ui
