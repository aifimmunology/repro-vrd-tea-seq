"""
Shared functions across different data apps.
"""
import argparse
from os import PathLike
from textwrap import wrap
from base64 import b64encode, b64decode
from io import BytesIO

from dash import dcc, Patch

import polars as pl

import plotly.graph_objects as go
import plotly.io as pio
from plotly.graph_objects import Figure

def setup_argparse():
    argparse_help = """
Allen Dash Modules (ADM) is a Python package that helps developers at the Allen Institute:
    1) write Dash apps more efficiently via pre-made classes, and 
    2) generalizes their Dash apps to multiple datasets using a onfigurable JSON file.
    """

    parser = argparse.ArgumentParser(
        description=argparse_help,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-c", "--config", required=False, help="Path to an ADM-compatible configuration file in JSON format.")
    parser.add_argument("-s", "--schema", required=False, help="Path to a JSON schema file to validate the above configuration file.")
    parser = parser.parse_args()

    return parser

def highlight_selected_gene(
        df: pl.DataFrame,
        fig: Figure,
        feature: str,
        gene: str,
        x_col: str,
        y_col: str,
        color_map: dict,
        title: str | None = None
    ) -> Figure:
    """
    Given a gene, create a new trace with a black border and text annotation above new point.
    Used for the following classes:
        * VolcanoPlot
        * MAPlot
    
    """

    # subset dataframe to selected gene
    gene_df=df.filter(pl.col(feature)==gene)

    if gene_df.shape[0] == 0:
        return fig

    # identify point color based on label
    gene_regulation_direction=color_map[gene_df["Category"][0]]

    # selected point information
    x=gene_df[x_col]
    y=gene_df[y_col]
    text=gene_df[feature][0]

    # add highlighted point
    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            marker_size=12,
            marker_color=gene_regulation_direction,
            marker_line_width=2,
            showlegend=False
        )
    )

    # add text above selected point
    fig.add_annotation(
        x=x[0],
        y=y[0],
        xref='x',
        yref='y',
        text=text,
        font=dict(
            size=20,
            weight="bold"
        ),
        yshift=20,
        showarrow=False
    )

    fig.update_layout(
        title=dict(text=f"Selected Gene: <b>{gene}</b>")
    )

    return fig

def highlight_gene_members(
        df: pl.DataFrame,
        fig: Figure,
        feature: str,
        genes: str,
        x_col: str,
        y_col: str,
        color_map: dict
    ) -> Figure:
    """
    Subset a dataframe given a list of genes, then highlight the points with a purple border.
    Used for the following classes:
        * VolcanoPlot
        * MAPlot

    """

    # filter df for gene members in a pathway. often a df with multiple rows.
    filtered_df = df.filter(pl.col(feature).is_in(genes))

    # return original figure if filtered df is empty
    if filtered_df.shape[0] == 0:
        return fig

    # define color of points
    marker_colors=[color_map[i] for i in filtered_df["Category"]]

    # modify figure
    fig.add_traces(
        go.Scattergl(
            x=filtered_df[x_col],
            y=filtered_df[y_col],
            marker_color=marker_colors,
            mode="markers",
            marker_size=12,
            marker_line_width=2,
            marker_line_color="#8b1ac4",
            marker_opacity=0.9,
            showlegend=False
        )
    )

    return fig

def highlight_point_with_label(
        df: pl.DataFrame,
        fig: Figure,
        feature_col: str,
        feature_value: str,
        x_col: str,
        y_col: str,
        color_map: dict,
        title = True,
        **kwargs
    ) -> Figure:
    """
    Annotate a point via an arrow and a text label that can wrap newlines if the label is too long.
    This is accomplished by providing a DF, subset by column feature_col that matches feature_value.
    Used for the following classes:
        * VolcanoPlot (GSEA)
        * VolcanoPlot (GEX by all conditions)

    """

    # prepare function call
    filtered_df = df.filter(pl.col(feature_col) == feature_value)

    # raise error message if empty df
    if filtered_df.shape[0] == 0:
        fig.layout.title = f"Results for <b>{feature_value}</b> are unavailable"
        return fig

    # skip action if no coordinates provided or any null values
    if filtered_df.select([x_col,  y_col]).drop_nulls().shape[0] == 0:
        fig.layout.title = f"Warning: no gene expression in <b>{feature_value}</b>"
        return fig

    # add highlighted point
    fig.add_trace(
        go.Scattergl(
            x=filtered_df[x_col],
            y=filtered_df[y_col],
            mode="markers",
            marker_size=12,
            marker_color=color_map[filtered_df["Category"][0]],
            marker_line_width=2,
            showlegend=False
        )
    )

    # arrow and text label
    feature_label=feature_value
    in_plot_label=feature_label
    yshift=20

    # wrap feature_label if too long
    if len(feature_value) > 40:
        in_plot_label = "<br>".join(wrap(feature_label, width=40))
        yshift=yshift*in_plot_label.count("<br>")

    fig.add_annotation(
        x=filtered_df[x_col][0],
        y=filtered_df[y_col][0],
        xref='x',
        yref='y',
        text=in_plot_label,
        font=dict(
            size=20,
            weight="bold"
        ),
        xshift=-1,
        yshift=7,

        showarrow=True,
        arrowhead=2,
        arrowwidth=2,
        ax=-10 if filtered_df[x_col][0] > 0 else 10,
        ay=-40,

        **kwargs
    )

    # wrap plot title with selected pathway
    if title is False:
        fig.update_layout(
            title=dict(text="")
        )
    else:
        if len(feature_label) > 50:
            title = "<br>".join(wrap(feature_label, width=50))

        fig.update_layout(
            title=dict(
                text=f"Selected Feature: <b>{feature_label}</b>")
            )

    return fig

def hex_to_rgba(hex_color, alpha=0.5) -> str:
    """
    Converts a hexadecimal color code to an RGBA string.

    Args:
        hex_color (str): The hexadecimal color code (e.g., '#FFFFFF', 'FFFFFF').
        alpha (float, optional): The alpha transparency value (0.0 to 1.0). Defaults to 0.5.

    Returns:
        str: The RGBA color string (e.g., 'rgba(255, 255, 255, 0.5)').
    """
    hex_color = hex_color.lstrip('#')  # Remove '#' if present

    # Ensure the hex color is 6 characters long (RRGGBB)
    if len(hex_color) == 3:
        hex_color = ''.join([c*2 for c in hex_color])
    elif len(hex_color) != 6:
        raise ValueError("Invalid hex color format. Must be 3 or 6 characters.")

    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)

    return f'rgba({r}, {g}, {b}, {alpha})'

def load_markdown(markdown_file: PathLike) -> dcc.Markdown:
    """
    This function reads the contents of a markdown file and returns a Dash Markdown component.

    Consumed by the following classes: HelpPanel, StartupModal

    Parameters
    ----------
    markdown_file : PathLike
        Path to the markdown file containing help content.

    Returns
    -------
    dash.dcc.Markdown
        A Markdown component containing the file's content with MathJax enabled.
    """

    # load markdown file
    with open(markdown_file, "r", encoding="utf-8") as file:
        page_content = file.read()

    # Convert to Markdown component with MathJax and className
    markdown_content = dcc.Markdown(
        page_content,
        mathjax=True,
        link_target="_blank",
        className='markdown-body'
    )

    return markdown_content

def make_ag_grid_column_filters(
        df: pl.LazyFrame | pl.DataFrame,
        exclude: list[str] | None = None
    ):

    # collect column names and types whether source comes
    # from a polars lazyframe or dataframe
    if isinstance(df, pl.LazyFrame):
        dtypes = df.collect_schema().dtypes()
        columns = df.collect_schema().names()
    elif isinstance(df, pl.DataFrame):
        dtypes = df.dtypes
        columns = df.columns
    else:
        exit("ERROR: input dataframe is not a polars lazyframe or dataframe")

    fields = []

    for dtype, col in zip(dtypes, columns):

        if exclude is not None:
            if col in exclude:
                continue

        # default column behavior
        column_configs = {
            'field': col,
            'sortable': True,
            'flex': 1,
            "filterParams": {
                "buttons": ["apply", "reset", "clear"], # Include "apply"
                "closeOnApply": True # Optional: closes the filter popup on apply
                }
            }

        # string columns should have filters for text
        if dtype == pl.String:
            column_configs['filter'] = 'agTextColumnFilter'
            fields.append(column_configs)

        # float columns should have filters for numbers
        elif dtype == pl.Float32:

            column_configs['filter'] = 'agNumberColumnFilter'

            # additionally, float columns with pvalues should use scientific notation with 3 sig figs,
            # unless the column is log or natural log transformed...
            # otherwise for any float column, display results with 4 sig figs
            column_name = col.lower()
            if "pval" in column_name or "adj" in column_name:
                if "log" in column_name or "ln" in column_name:
                    column_configs["valueFormatter"] = {
                        "function": """d3.format(",.4f")(params.value)"""
                    }
                else:
                    column_configs["valueFormatter"] = {
                        "function": """d3.format(",.4e")(params.value)"""
                    }
            else:
                column_configs["valueFormatter"] = {
                    "function": """d3.format(",.4f")(params.value)"""
                }

            fields.append(column_configs)

    return fields

def serialize_lazyframe(lazy_frame: pl.LazyFrame) -> str:
    """
    Given a lazyframe, turn it into a serializable string to be passed between callbacks.
    """
    return b64encode(lazy_frame.serialize()).decode('utf-8')

def deserialize_lazyframe(string: str) -> pl.LazyFrame:
    """
    Given a string, turn it into a lazyframe to be passed between callbacks.
    """
    return pl.LazyFrame.deserialize(BytesIO(b64decode(string)))

def patch_figure_theme(light: bool) -> Patch:

    template = "plotly_white" if light else "plotly_dark"

    fig = Patch()
    fig['layout']['template'] = pio.templates[template]
    return fig
