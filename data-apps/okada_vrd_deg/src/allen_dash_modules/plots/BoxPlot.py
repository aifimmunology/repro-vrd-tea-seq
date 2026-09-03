from __future__ import annotations
import sys

from dash_bootstrap_components import Card, CardBody, CardHeader
from dash.dcc import Graph, Loading

import polars as pl

import plotly.graph_objects as go
from itertools import cycle

def assign_shapes_to_samples(samples, shape_codes = [0, 1, 2, 3, 4, 13, 14, 17, 18, 21, 22]):
    """
    Assigns a shape code to each sample in a cyclic manner.
    
    Parameters:
        samples (list): A list of sample names/IDs.
        shape_codes (list): A list of Plotly marker symbol codes.
    
    Returns:
        dict: Mapping of sample -> shape_code.
    """
    shape_cycle = cycle(shape_codes)
    return {sample: next(shape_cycle) for sample in samples}

class BoxPlot:
    def __init__(self, id, title):
        self.id = id
        self.build_ui(title=title)

    def plot(
        self,
        df:                         pl.DataFrame,
        title:                      str,
        template:                   str,
        condition:                  str,
        perturbation:               bool        = False,
        perturbation_subject_label: str | None  = None,
        perturbation_control_label: str | None  = None
    ):
        """
        Example input data...
        ┌───────┬─────────┬──────────┬───────────────────────┬────────────┐
        │ Gene  ┆ Donor   ┆ Cytokine ┆ Cell Type             ┆ Expression │
        │ ---   ┆ ---     ┆ ---      ┆ ---                   ┆ ---        │
        │ str   ┆ str     ┆ str      ┆ str                   ┆ f64        │
        ╞═══════╪═════════╪══════════╪═══════════════════════╪════════════╡
        │ PDE4B ┆ Donor1  ┆ 4-1BBL   ┆ B-Intermediate/Memory ┆ 1.241975   │
        │ PDE4B ┆ Donor1  ┆ PBS      ┆ B-Intermediate/Memory ┆ 1.214031   │
        │ PDE4B ┆ Donor10 ┆ 4-1BBL   ┆ B-Intermediate/Memory ┆ 1.015893   │
        │ PDE4B ┆ Donor10 ┆ PBS      ┆ B-Intermediate/Memory ┆ 0.855644   │
        │ PDE4B ┆ Donor11 ┆ 4-1BBL   ┆ B-Intermediate/Memory ┆ 0.788048   │
        │ PDE4B ┆ Donor11 ┆ PBS      ┆ B-Intermediate/Memory ┆ 1.020474   │
        └───────┴─────────┴──────────┴───────────────────────┴────────────┘

        Gene and other independent variables will be a constant.

        To make a paired boxplot, make sure `paired` is True and `condition` is set to the column ("Donor") that repeats twice.
        Then designate one independent variable as the condition/treatment group (Cytokine in this case).

        If this is a perturbation experiment, set perturbation=True and define the perturbation_control_label='PBS'.
        This will draw lines and markers that connect across box plots.

        NOTE: only 2-condition perturbations are supported at the moment.
        """

        # instantiate constants
        fig = go.Figure()

        if not perturbation:

            # 1. Create the boxplot (empty fill with black points).
            for condition_group in df[condition].unique():
                fig.add_trace(
                    go.Box(
                        y=df.filter(pl.col(condition) == condition_group)["Expression"].to_list(),
                        x=[condition_group] * df.filter(pl.col(condition) == condition_group).height,
                        name=condition_group,
                        boxpoints="all",
                        pointpos=0,
                        fillcolor='rgba(0,0,0,0)',
                        line_color='#808080',
                        line_width=2,
                        showlegend=False
                    )
                )

        # 2. If perturbation study, add paired points via line and marker.
        else:

            # generate similar boxplot as before without points
            for condition_group in df[condition].unique():
                fig.add_trace(
                    go.Box(
                        y=df.filter(pl.col(condition) == condition_group)["Expression"].to_list(),
                        x=[condition_group] * df.filter(pl.col(condition) == condition_group).height,
                        name=condition_group,
                        boxpoints=None,
                        fillcolor='rgba(0,0,0,0)',
                        line_color='#C7C7C7',
                        line_width=2,
                        showlegend=False
                    )
                )

            if perturbation_control_label in df[condition].unique().to_list():

                treatment=condition
                perturbed_condition_label = [ i for i in df[treatment].unique().to_list() if i != perturbation_control_label ]
                order_of_labels = [perturbation_control_label, perturbed_condition_label]

                if perturbation_subject_label is not None:
                      paired_subjects=df[perturbation_subject_label].unique().sort()
                else:
                    print(df)
                    sys.exit(f"ERROR: perturbation_subject_label {perturbation_subject_label} not found in dataframe")

                # color_list = px.colors.qualitative.Set3
                color_list = ['#0d6efd', '#d63384', '#fd7e14', '#ffc107', '#198754', '#20c997', '#d63384', '#5286b0', '#A86F73', '#563210', '#90CBB0', '#FE0000']
                color_map = dict(zip(paired_subjects, color_list * ((len(paired_subjects) // 10) + 1)))
                symbol_map = assign_shapes_to_samples(paired_subjects)

                for subject in paired_subjects:
                    subject_df = df.filter(pl.col(perturbation_subject_label) == subject)
                    if subject_df.height == 2:
                        treatment_array = subject_df[condition].to_list() # categorical list
                        expressions_array = subject_df["Expression"].to_list() # numerical list
                        fig.add_trace(
                            go.Scatter(
                                x=treatment_array,
                                y=expressions_array,
                                mode="lines+markers",
                                name=subject,
                                line=dict(color="#C7C7C7", width=2),
                                marker=dict(size=16, color=color_map[subject]),
                                marker_symbol=symbol_map[subject],
                                hoverinfo="text",
                                text=[f"{subject} ({c})" for c in treatment_array],
                                showlegend=True
                            )
                        )
                fig.update_xaxes(categoryarray=order_of_labels)
            else:
                print(df)
                sys.exit(f"ERROR: {condition} was not found in the condition column.")

        fig.update_xaxes(
            showline = True,
            linecolor = 'black',
            gridcolor = 'lightgrey',
            automargin=True
        )

        fig.update_yaxes(
            automargin=True
        )

        fig.update_layout(

            title=title,
            xaxis_title="Treatment",
            yaxis_title="Expression",

            hoverlabel=dict(bordercolor="rgb(0,0,0,0)"),
            clickmode = 'event',
            font_family='Roboto Condensed',
            font_size=14,

            autosize = True,

            xaxis=dict(title_standoff=10),
            yaxis=dict(title_standoff=25),

            template = template
        )

        return fig

    def build_ui(self, title):

        ui = Card([
            CardHeader(title),
            CardBody([
                Loading(
                    Graph(
                        id = self.id,
                        figure = {},
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
                    "filter": "blur(2px)"},
                type='circle',
                delay_show=2000
                )  # end Loading
            ]) # end CardBody
        ],
        style={
            'height': '100%',
            'width': '100%'
        } # end Card children
    ) # end Card

        self.ui = ui
