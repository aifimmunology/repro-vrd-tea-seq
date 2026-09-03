app_version = '2.0.0'  #20230907

from dash import Dash, dcc, html, Input, Output, State, MATCH, clientside_callback
from dash_iconify import DashIconify

import src.allen_dash_modules as adm

# Set Up
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
import plotly.subplots
import dash_daq as daq # import needed, for numeric input
import plotly

import polars as pl
import json

import math
import numpy as np
import pandas as pd
import re
from sklearn.cluster import AgglomerativeClustering

# Gene info
import mygene  # gene info
from dash import dash_table  # gene info

import logging
from warnings import warn

### Official HISE colors
colors = {
    "--bs-blue": "#0d6efd",
    "--bs-indigo": "#6610f2",
    "--bs-purple": "#6f42c1",
    "--bs-pink": "#d63384",
    "--bs-red": "#dc3545",
    "--bs-orange": "#fd7e14",
    "--bs-yellow": "#ffc107",
    "--bs-green": "#198754",
    "--bs-teal": "#20c997",
    "--bs-cyan": "#0dcaf0",
    "--bs-white": "#FFFFFF",
    "--bs-gray": "#6c757d",
    "--bs-gray-dark": "#343a40",
    "--bs-gray-100": "#f8f9fa",
    "--bs-gray-200": "#e9ecef",
    "--bs-gray-300": "#dee2e6",
    "--bs-gray-400": "#ced4da",
    "--bs-gray-500": "#adb5bd",
    "--bs-gray-600": "#6c757d",
    "--bs-gray-700": "#495057",
    "--bs-gray-800": "#343a40",
    "--bs-gray-900": "#212529",
    "--bs-hise-blue-1": "#003056",
    "--bs-hise-blue-2": "#325876",
    "--bs-hise-blue-3": "#5286b0",
    "--bs-hise-blue-4": "#71899c",
    "--bs-hise-blue-5": "#b4c3cf",
    "--bs-hise-blue-6": "#d9e0e6",
    "--bs-hise-teal-1": "#33B0C8",
    "--bs-hise-teal-2": "#76CFE0",
    "--bs-hise-teal-3": "#DEF2F6",
    "--bs-hise-grey-1": "#272D3B",
    "--bs-hise-grey-2": "#3E3E3E",
    "--bs-hise-grey-3": "#616161",
    "--bs-hise-grey-4": "#707070",
    "--bs-hise-grey-5": "#ECEDEE",
    "--bs-hise-grey-6": "#FBFBFB",
    "--bs-hise-green-1": "#E3EADA",
    "--bs-hise-green-2": "#A0C572",
    "--bs-hise-green-3": "#94BC62",
    "--bs-hise-green-4": "#4AD991",
    "--bs-aifi-new-001": "#003057",
    "--bs-aifi-new-002": "#5da7e5",
    "--bs-aifi-new-003": "#74A03E",
    "--bs-aifi-new-004": "#f4a261",
    "--bs-aifi-new-005": "#e76f51",
    "--bs-aifi-new-006": "#FFFFD0",
}

# create header object
header = adm.AppHeader(
    toggle_help=True,
    toggle_settings=False,
    toggle_theme=False
)

# build app header
header.ui = header.build_app_header(
    title="T Cell RVd Differential Expression Explorer",
    institute_url="https://alleninstitute.org/division/immunology/",
    study_url="https://apps.allenimmunology.org/aifi/insights/tcell-vrd/"
)

# build help page by providing a markdown file
header.help_page.ui = header.help_page.build_help_panel(
    markdown_file="data/app_info.md"
)

# establish stylesheets
stylesheets = [
    "https://cdn.jsdelivr.net/npm/bootstrap@5.3.6/dist/css/bootstrap.min.css",
]

# build dash app
app = Dash(external_stylesheets=stylesheets, update_title="")

# swap favicon to match user's light/dark theme preference.
# dark theme => white svg, light theme => black svg.
clientside_callback(
    """
    (_) => {
        const setFavicon = (isDark) => {
            const href = isDark
                ? "/assets/images/AI_lens_white.svg"
                : "/assets/images/AI_lens_black.svg";
            let link = document.querySelector("link[rel~='icon']");
            if (!link) {
                link = document.createElement("link");
                link.rel = "icon";
                document.head.appendChild(link);
            }
            link.type = "image/svg+xml";
            link.href = href;
        };
        const mql = window.matchMedia('(prefers-color-scheme: dark)');
        setFavicon(mql.matches);
        if (!window.__faviconThemeListener) {
            window.__faviconThemeListener = true;
            mql.addEventListener("change", (e) => setFavicon(e.matches));
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("ui-header", "id"),
    Input("ui-header", "id"),
)


## Functions

### volcano

def assign_change(p_val, p_cutoff, es, es_cutoff, coding = ['dn','nc','up']):
    change = np.repeat(coding[1], len(p_val), axis=0)
    iUp =  [i for i, (l1, l2) in enumerate(zip(es, p_val)) if l1 is not None and l2 is not None and l1 >= es_cutoff and l2 <= p_cutoff]
    if len(iUp) > 0 :  
        change[iUp] = coding[2]
    iDown = [i for i, (l1, l2) in enumerate(zip(es, p_val)) if l1 is not None and l2 is not None and l1 <= -1*es_cutoff and l2 <= p_cutoff]
    if len(iDown) > 0 : 
        change[iDown] = coding[0]
        
    return(change)
    
# is there a faster way to do this using just go? 
def build_volcano(
    de_df, 
    colorcol = None, 
    p_name = 'adjP',
    es_name = 'logFC', 
    alpha = 0.05,
    es_cutoff = 0.2, 
    hovername = 'primerid',
    color_dictionary = {'up' : '#e3480b','dn' : '#0b7be3','nc' : '#808080'},
    opacity_val = 0.3,
    highlight_opacity = 1,
    highlight_color = '#000000',
    logtransp = True,
    min_p = 1e-200,
    height = 500,
    width = 500,
    title = '',
    title_size = 20,
    highlight = None,
    size_val = 5,
    highlight_size = 7
    ):
    
    if colorcol is None:  
        colorcol = 'Change'
        de_df.loc[:,'Change'] = assign_change(
            p_val = de_df.loc[:,p_name].tolist(), 
            p_cutoff= alpha, 
            es = de_df.loc[:,es_name].tolist(), 
            es_cutoff = es_cutoff
        )
        
    if logtransp:
        de_df.loc[:,p_name] = [x if x >= min_p else min_p for x in de_df.loc[:,p_name]]  # ensure no 0 values before log transform
        de_df.loc[:,'-1*log10('+p_name+')'] =  pd.DataFrame.apply(de_df, func = lambda row: -1*math.log10(row[p_name]), axis=1)
        yname = '-1*log10('+p_name+')'
    else:
        yname = p_name
    
    if highlight is not None and highlight != 'All' and len(highlight) > 0:
        if type(highlight) != list:
            highlight = [highlight]
        
        if any(de_df.loc[:, hovername].isin(highlight)):
            for nm in highlight:
                if nm not in de_df.loc[:, hovername].values:
                    print('Warning: Feature(s) "' + nm + '" not present in expected column ' + hovername + '.')
            hlt_index = de_df[hovername].isin(highlight)
            hlt_subset1 = de_df.loc[hlt_index,[es_name, yname,hovername]]
            hlt_subset2 = de_df.loc[hlt_index,]
            de_df = de_df[~de_df[hovername].isin(highlight)]
        else:
            highlight_str = highlight[0]
            if len(highlight) > 1:
                for nm in highlight[1:]:
                    highlight_str = highlight_str + ', ' + nm
            print('Warning: Feature(s) ' + highlight_str + ' not present in expected column ' + hovername + '. Showing all features.')
    
    fig = px.scatter(
        de_df, 
        x=es_name, 
        y=yname,
        hover_name = hovername,
        hover_data = [es_name, yname],
        color = colorcol,
        color_discrete_map = color_dictionary,
        opacity = opacity_val,
        title = title,
        render_mode='webgl'
    ).update_traces(
        marker_size=size_val
    )
    
    if 'hlt_subset1' in locals():
        for x,y,z in hlt_subset1.itertuples(index=False):
            fig.add_annotation(
                x=x, 
                y=y,
                text=z,
                showarrow=False,
                yshift=10,
            )

            fig.add_traces(     
                px.scatter(
                    hlt_subset2,
                    x=es_name, 
                    y=yname,
                    hover_data = [colorcol, es_name, yname],
                    hover_name = hovername,
                    render_mode='webgl'  # force webgl so the fig.data order is preserved, with highlighted on top (last)
                ).update_traces(
                    marker_size=highlight_size,
                    marker_color=highlight_color, 
                    marker_opacity = 1,
                    name = 'highlighted'
                ).data
            )
    
    fig.update_layout(
        autosize=False,
        width=width,
        height=height,
        clickmode='event+select',
        title={'font':{'size':title_size}}
    )

    return fig

# ### Gene Info 

def _format_field(value, field, qry_res):
    value_dict = get_output_field(qry_res, field)

    if value is None:
        out_text = 'No gene selected'
    elif value in value_dict.keys():
        out_text = value_dict[value]
        if type(out_text) == list:
            out_text = '; '.join(out_text)
        if out_text is None:
            out_text = 'No Result'
    else:
        out_text = 'No Result'
    return out_text


def make_go_table(value, index, qry_res = None):
    go_col_list = ['evidence', 'gocategory','id', 'pubmed','qualifier', 'term']
    go_types_dict = {
        'evidence': 'text',
        'gocategory': 'text',
        'id': 'text',
        'pubmed': 'text',
        'qualifier': 'text',
        'term': 'text'
    }
    
    if qry_res is None:
        qry_res =  query_genes(value, 'go')
    
    go_dict = get_output_field(qry_res, 'go') 

    if (value is None) or (value not in go_dict.keys()) or (go_dict[value] is None):
        table_data = None
        n_records = 0
        countstring = 'No matching records'
    else:
        table_data = compile_go_df(go_dict[value])
        table_data = table_data.to_dict('records')
        n_records = len(table_data)
        countstring = 'Total records: ' + str(n_records) #+ '/Total pages: ' + str(math.ceil(n_records/5))
    
    go_div = html.Div(
        id={
            'type':'gene-summary-container2',
            'index':index,
        },
        className='col-6',
        children = [
            html.Br(),
            html.H3('Gene GO References'),
            html.Div(id = 'go_table_recordcount', children = countstring),
            html.Div('Select records per page'),
            dcc.Dropdown(    
            id={
                'type':'select_page_size',
                'index':index
            }, 
            options=[{'label': '5', 'value': 5}, {'label': '10', 'value': 10}, {'label': '25', 'value': 25}],    
            value=5    
            ),
            html.Div(
                id={
                    'type':'go_table_recordcount',
                    'index':index,
                },
                children = []
            ),
            html.Div(
                id={
                    'type':'go_scroll_container',
                    'index':index,
                },
                style = {'overflow-x': 'auto','overflow-y': 'auto'},
                children = [
                dash_table.DataTable(
                    id={
                        'type':'gene-go-container',
                        'index':index
                    },
                    style_cell = {'font_family':'system-ui'},
                    data = table_data,
                    columns=[{'id': x, 'name': x, 'type': go_types_dict[x]} for x in go_col_list],
                    page_current=0,
                    page_size = 5,
                    page_action='native',
                    filter_action='native',
                    filter_query='',
                    sort_action='native',
                    sort_mode='multi',
                    sort_by=[])
                ]
            )  # end row 1 column 2 scrolling container
        ]  # end row 1 column 2 children
    )  # end row 1 column 2
    
    return go_div

def query_genes(value, myfields=['name','alias','generif','genomic_pos','go','other_names','pathway','summary']):
    mg = mygene.MyGeneInfo()
    myfields = ['name','alias','generif','genomic_pos','go','other_names','pathway','summary']
    qry_res =  mg.querymany(value, scopes='symbol', fields = myfields, species='human')
    
    return(qry_res)


def make_geneinfo_container(value, index, qry_res = None, myfields=['name','alias','generif','genomic_pos','go','other_names','pathway','summary']):
    if qry_res is None:
        qry_res = query_genes(value, myfields)
    
    gene_div = html.Div(
        className = 'row',
        style = {'background-color': 'white', 'padding': '10px'},
        children = [
            html.H2('Selected Gene: {}'.format(value)),
            html.Div(
                id={
                    'type':'gene-summary-container',
                    'index':index,
                },
                className='col-6',
                children = [
                    html.Br(),
                    html.H3('Gene Name'),
                    html.Div(
                        children = _format_field(value, 'name', qry_res)
                    ),
                    html.Br(),
                    html.H3('Other Names'),
                    html.Div(
                        children = _format_field(value, 'other_names', qry_res)
                    ),
                    html.Br(),
                    html.H3('Gene Aliases'),
                    html.Div(
                        children = _format_field(value, 'alias', qry_res)
                    ),
                    html.Br(),
                    html.H3('Gene Summary'),
                    html.Div(
                        children = _format_field(value, 'summary', qry_res)
                    )
                ]
            ),
            make_go_table(value, index, qry_res)
        ] # end row 1 children
    )
    
    return gene_div

# Get value of field in mygene query output. 
#' Expecting first level to be list query results, second level to be 
#' result for specific field for that query (gene)
# ls = list of query results from mygene. for example, output of mygene.querymany() where returnall=False, or
#   mygene.queryman()['out'] if returnall=True
#' field = which field from returned query to return
#' namefield = which field in output is the key, default 'query' (genesymbol)
def get_output_field(ls, field, namefield = 'query'):
    out_dict = {}
    missing_count = 0
    for i in range(len(ls)):
        if(namefield is None): # ie if you queried genes using query() and there is no query field
            if field in ls[i].keys():
                val = ls[i][field]
                out_dict[i-missing_count] = val
            else:
                missing_count += 1
        else: # ie if you queried genes using querymany()
            nm = ls[i][namefield]
            if field in ls[i].keys():
                val = ls[i][field]
                if nm in out_dict.keys():  # there were multiple hits for the same query, merge output
                    if out_dict[nm] is None:
                        out_dict[nm] = val
                    if type(out_dict[nm]) == 'string':
                        out_dict[nm] = out_dict[nm] + '\n' + val
                    elif type(out_dict[nm]) == 'list':
                        out_dict[nm] = out_dict[nm].extend(val)
                else:
                    out_dict[nm] = val
            else:
                out_dict[nm] = None
    return(out_dict)

# Compile the gene ontology mygene output into a data frame
def compile_go_df(go_dict):
    out_df_ls = []
    
    for key in go_dict:
        
        sub_d = go_dict[key]
            
        if not type(sub_d) == list:
            sub_d = [sub_d]  # handle case if one go type only had one entry (returns single dictionary instead of list of dictionaries

        df = pd.DataFrame.from_dict(sub_d, orient = 'columns')
        if 'category' in df.columns :
            df.rename(columns={'category' : 'gocategory'}, inplace = True)  # GO MF result has 'category' column instead of 'gocategory', make consistent
        
        if (len(out_df_ls) < 1):
            out_df_ls = [df]
        else:
            out_df_ls.append(df) 
    
    dfout = pd.concat(out_df_ls)
    for colname in dfout.columns:
        dfout[colname]= dfout[colname].map(str)
    
    return  dfout

def add_to_list(ls, new_ls):
    for x in new_ls:
        if x not in ls:
            ls.extend([x])
    
    return(ls)


### Subplot Metadata annotations
# take a matplotlib figure and return a string encoding that can be rendered in dash using html.Img()
def fig_to_encoding(fig, save_fmt= 'jpg', out_fmt='png', decoding = 'ascii'):
    my_stringIObytes = io.BytesIO()
    plt.savefig(my_stringIObytes, format=save_fmt)
    my_stringIObytes.seek(0)
    encoded = base64.b64encode(my_stringIObytes.read()).decode(decoding).replace('\n', '')
    out_string = 'data:image/{};base64,{}'.format(out_fmt, encoded)
    my_stringIObytes.close()
    
    return(out_string)

def discrete_colorscale(bvals, colors):
    """Generates a interval-based discrete colorscale based on input bin values and colors
    
    Parameters
    ----------
    bvals : list
        Numeric values bounding intervals/ranges of interest.
    colors : list
        Colorcodes (mapping to rgb, rgba, hex, hsl, hsv, or named color string) for bins defined by bvals. Should be of length len(bvals)-1.
    
    Returns
    -------
    list
        The plotly discrete colorscale that can be passed to a plotly graph object as colorscale parameter.
    
    Notes
    -----
    source: https://chart-studio.plotly.com/~empet/15229.embed
    
    Examples
    --------
    >>> discrete_colorscale([0,1,5,9,10], ['red','orange','yellow','green'])
    [
        [0.0, 'red'],
        [0.1, 'red'],
        [0.1, 'orange'],
        [0.5, 'orange'],
        [0.5, 'yellow'],
        [0.9, 'yellow'],
        [0.9, 'green'],
        [1.0, 'green']
    ]
    
    """
    if len(bvals) != len(colors)+1:
        raise ValueError('len(boundary values) should be equal to  len(colors)+1')
    bvals = sorted(bvals)     
    nvals = [(v-bvals[0])/(bvals[-1]-bvals[0]) for v in bvals]  #normalized values
    dcolorscale = [] #discrete colorscale
    for k in range(len(colors)):
        dcolorscale.extend([[nvals[k], colors[k]], [nvals[k+1], colors[k]]])
    
    return dcolorscale    

def make_cat_colorscale(meta, var, color_dict):
    """Generates a interval-based discrete colorscale based on a dictionary of categories and colors.
    
    Uses the dictionary to generate color boundaries based on category levels.
    
    Parameters
    ----------
    meta : pandas.DataFrame
        Data frame of metadata containing the variable of interest. Used for determining variable category levels.
    var : str 
        Metadata variable name
    color_dict: dictionary
        Key are levels of variable var, and values are colors (Ie 'red', 'rgb(13, 8, 135)', etc)
    
    Returns
    -------
    list
        A plotly discrete colorscale that can be passed to a plotly graph object as colorscale parameter.
    
    Examples
    --------
    >>> test_meta = pd.DataFrame(
        [
            [2,'two'],
            [1,'one'],
            [3,'one'],
            [4,'four'],
            [5,'three']
        ], columns = ['Item','Group'])
    >>> test_meta['Group'] = pd.Categorical(test_meta['Group'], categories =['one','two','three','four'])
    >>>  make_cat_colorscale(test_meta, 'Group', dict(one = 'red', two = 'green', three = 'orange', four = 'purple'))
    [
        [0.0, 'red'],
        [0.25, 'red'],
        [0.25, 'green'],
        [0.5, 'green'],
        [0.5, 'orange'],
        [0.75, 'orange'],
        [0.75, 'purple'],
        [1.0, 'purple']
    ]
    """
    cat_lev = meta[var].cat.categories
    bvals = range(len(cat_lev))
    bvals=np.append(bvals, [len(bvals)])
    cols = [color_dict[x] for x in cat_lev]
    colorscale = discrete_colorscale(bvals, cols)
    return(colorscale)

def _calc_subplot_row_heights(nfeat, nmeta):
    """Calculate rows of plotly subplots for annotated heatmap
    
    Parameters
    ----------
    nfeat : int
        Number of features plotted on the main heatmap
    nmeta : int
        Number of metadata annotations plotted (as separate rows)s

    Returns
    -------
    list 
        list of heights for each plot, in inches
    
    Examples
    --------
    >>> _calc_subplot_row_heights(100, 3)
    [
        0.009671179883945842,
        0.009671179883945842,
        0.009671179883945842,
        0.9671179883945842
    ]
    
    """
    tot_row = nfeat + 0.1 + nmeta*1.1
    rowval = 1/tot_row
    heights = pd.Series([rowval]).repeat(nmeta).tolist()
    heights.extend([nfeat*rowval])

    return(heights)

def _calc_legend_height(meta, meta_plot, height_per_level = 0.04, continuous_height = 0.2, matrix_val_height = 0.2, legend_font_size=12):
    """Calculate legend heights of plotly subplots for annotated heatmap
    
    Parameters
    ----------
    meta : pandas data frame
        Metadata table, one row per observation (ie per cell)
    meta_plot : list
        List of metadata columns in meta to include in heatmap annotations
    height_per_level: float, default 0.04
        Colorbar height allocated to each level of a categorical variable in the legend. Fractional value between 0 and 1, but should be <= 1/nlevels.
    continuous_height : float, default 0.2
        Entire colorbar height for a continuous variable in the legend. Fractional value between 0 and 1, but should be <= 1.    
    matrix_val_height : float, default 0.2
        Entire colorbar height forthe main matrix measurement in the legend. Fractional value between 0 and 1, but should be <= 1.        
    legend_font_size : int, default 12
    
    Returns
    -------
    list 
        list of heights for each legend colorbar, in fractions. The last item will be the colorbar for the main matrix heatmap.
    
    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> test_meta = pd.DataFrame(
        [
            [2,'two','banana',8.0],
            [1,'one','apple',5.9],
            [3,'one','banana',2.3],
            [4,'four','apple',12.4],
            [5,'three','apple',2.1]
        ], 
        columns = ['Item','Group','Fruit','Weight']
    )
    >>> test_meta['Group'] = pd.Categorical(test_meta['Group'], categories =['one','two','three','four']) 
    >>> test_meta['Fruit'] = pd.Categorical(test_meta['Fruit']) 
    >>> _calc_legend_height(test_meta, ['Group','Fruit','Weight'])
    [0.32, 0.24, 0.2, 0.2]
    
    """
    heights = []
    for var in meta_plot:
        if isinstance(meta[var].dtype, pd.CategoricalDtype):
            levels = meta[var].cat.categories.values
            nlevels = len(levels)
            height = height_per_level*(nlevels+4*legend_font_size/12)
            heights.extend([height])
        else:
            heights.extend([continuous_height])
            
    # Add matrix legend at end if included
    if matrix_val_height is not None:
        heights.extend([matrix_val_height])
    return(heights)


def _calc_meta_legend_locs(meta, meta_plot, height_per_level = 0.03, continuous_height = 0.2, matrix_val_height=0.2, xshift = 0.1, yshift = 0.01, legend_font_size =12):
    """Calculate the legend positions for each heatmap subplot
    
    The legend colorbar placement will start at position 1, 1 at the upper right of the plotting field. Each subsequent colorbar will be placed yshift distance below
    the end of the previous colorbar unless it would fall out of frame, in which case it would be shifted right by xshift and placed at the top of the plot.
    
    Parameters
    ----------
    meta : pandas data frame
        Metadata table, one row per observation (ie per cell)
    meta_plot : list
        List of metadata columns in meta to include in heatmap annotations
    height_per_level: float, default 0.04
        Colorbar height allocated to each level of a categorical variable in the legend. Fractional value between 0 and 1, but should be <= 1/nlevels.
    continuous_height : float, default 0.2
        Entire colorbar height for a continuous variable in the legend. Fractional value between 0 and 1, but should be <= 1.    
    matrix_val_height : float, default 0.2
        Entire colorbar height forthe main matrix measurement in the legend. Fractional value between 0 and 1, but should be <= 1.   
    xshift : float, default 0.1
        Fractional value, horizontal spacing between starting locations of legend colorbar columns
    yshift : float, default 0.1
        Fractional value, vertical spacing between adjacent colorbars in legend
    legend_font_size : int, default 12
    
    Returns
    -------
    list 
        List of length 2. First item is a list of x positions and second item is a list of y positions corresponding to each legend colorbar.
    
    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> test_meta = pd.DataFrame(
        [
            [2,'two','banana',8.0],
            [1,'one','apple',5.9],
            [3,'one','banana',2.3],
            [4,'four','apple',12.4],
            [5,'three','apple',2.1]
        ], 
        columns = ['Item','Group','Fruit','Weight']
    )
    >>> test_meta['Group'] = pd.Categorical(test_meta['Group'], categories =['one','two','three','four']) 
    >>> test_meta['Fruit'] = pd.Categorical(test_meta['Fruit']) 
    >>> _calc_meta_legend_locs(test_meta, ['Group','Fruit','Weight'])
    [[1, 1, 1, 1], [1, 0.75, 0.56, 0.35000000000000003]]
    
    """
    legend_heights = _calc_legend_height(meta, meta_plot, height_per_level, continuous_height, matrix_val_height,legend_font_size)
    legend_y = [1]
    legend_x = [1]
    for i in range(len(legend_heights)):
        if i == 0:
            next
        else:
            yval = legend_y[i-1] - legend_heights[i-1]- yshift 
            if (yval - legend_heights[i]) < 0:
                # move to next column
                legend_y.extend([1])
                legend_x.extend([legend_x[i-1] + xshift])
            else:
                legend_y.extend([yval])
                legend_x.extend([legend_x[i-1]])
    
    return([legend_x, legend_y])

# for adding custom field(s) for hoverdata string (hover_template). must add customdata param to plot call. see format_heatmap_anno()
def _get_custom(data):
    out_str = ''
    for i in range(data.shape[1]):
        cname = data.columns[i]
        if pd.api.types.is_float_dtype(data[cname]):
            fmtstring = ':.3f'
        else:
            fmtstring = ''
        istring = cname + ': %{customdata['+ str(i) +']'+fmtstring+'}' + '<br>' 
        out_str = out_str + istring
    
    return(out_str)

def format_heatmap_anno(meta, var, xcol='barcodes', yval=None, colorscale='plasma', legendgrp='1', xloc=1, yloc=1, height_per_level=0.03, continuous_height=0.2, legend_font_size=12):
    """Creates a heatmap of a metadata feature for annotating alongside an expression heatmap
    
    Parameters
    ----------
    meta : pandas data frame
        Metadata table, one row per observation (ie per cell)
    var : str
        Metadata variable name to use on the z-axis (color values). If categorical, will plot discrete levels.
    xcol : str
        metadata variable name to use on the x-axis. Default 'barcodes' to plot each cell
    yval : str
        metadata variable name to use on the y-axis. Default None will use the 'var' column name as a single y-axis label.
    colorscale : str or dictionary, default='plasma'
        Either a custom dictionary of colors for each level of `var`, or a plotly-supported color palette. See
        `plotly.express.colors.named_colorscales()` for available palettes.
    legendgrp : str, default '1'
        Name of legend group the figure will belong to. For different subplots in a larger figure, use different group names. 
    xloc : float, optional
        X position value (between 0 and 1) of the plot legend in relation to the entire figure (not subplot). Default 1 for right side of plot.
    yloc : float, optional
        Y position value (between 0 and 1) of the plot legend in relation to the entire figure (not subplot). Default 1 for right side of plot.
    height_per_level : int, optional
        A fraction value, used for categorical variables only. The height alotted to each level of the variable in the legend colorbar. Should not exceed 1/nlevels. Default 0.03.
    continuous_height : float, optional
        A fraction value, for continuous variable only. The total height alotted to the colorbar in the legend. Should not exceed 1.
    legend_font_size : int, optional
        Fontsize for the legend text, default 12.
    
    Returns
    -------
    plotly.graph_objects.Heatmap object
        A heatmap object with one row and columns for each observation in meta. Colored by selected variable. Intended to be incorporated in a figure as a plotly subplot
        
    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> import plotly
    >>> test_meta = pd.DataFrame(
            [
                [2,'two','banana'],
                [1,'one','apple'],
                [3,'one','banana'],
                [4,'four','apple'],
                [5,'three','apple']
            ], 
            columns = ['Item','Group','Fruit']
        )
    >>> test_meta['Group'] = pd.Categorical(test_meta['Group'], categories =['one','two','three','four']) 
    >>> test_meta['Fruit'] = pd.Categorical(test_meta['Fruit']) 
    >>> group_cols = make_cat_colorscale(test_meta, 'Group', dict(one = 'red', two = 'green', three = 'orange', four = 'purple'))
    >>> fig = plotly.subplots.make_subplots(rows=2, cols=1, vertical_spacing = 0.01,shared_xaxes=True)   
    >>> hm1 = format_heatmap_anno(test_meta, 'Group', xcol='Item', yval=None, colorscale=group_cols, legendgrp='1',height_per_level=0.05, xloc=1, yloc=1, legend_font_size=12)
    >>> hm2 = format_heatmap_anno(test_meta, 'Fruit', xcol='Item', yval=None, colorscale='viridis', legendgrp='2', height_per_level=0.05, xloc=1, yloc=0.6, legend_font_size=12)
    >>> fig.append_trace(hm1, row = 1, col = 1)
    >>> fig.append_trace(hm2, row = 2, col = 1)
    >>> fig.update_layout(
        height=400,
        width=600)
    >>> fig.show()
    
    """
    
    varcol = var
    
    # y value, use variable name as a label
    if yval is None:
        yval = var
    yvals = pd.Series([yval]).repeat(meta.shape[0])
    
    metatype = meta[var].dtype.name
    if metatype == "category":
        # Make categorical column to fill by
        varcol = "__" + var + "_cat"
        meta[varcol] = meta[var].cat.codes
        
        bvals=np.unique(np.sort(meta[varcol]))
        bvals=np.append(bvals, [len(bvals)]) # add end to top bin
        if type(colorscale) == str:
            cols = px.colors.sample_colorscale(colorscale, len(bvals)-1)
            colorscale = discrete_colorscale(bvals, cols)
        elif type(colorscale) == dict:
            colorscale = make_cat_colorscale(meta, var, color_dict = colorscale)
            
        bvals = np.array(bvals)
        tickvals = [np.mean(bvals[k:k+2]) for k in range(len(bvals)-1)] #position with respect to bvals where ticktext is displayed
        labels = meta[var].cat.categories
        tmode = "array"
        cbar = dict(
            tickvals=tickvals, 
            ticktext=labels,
            len = height_per_level*(len(bvals)+3*legend_font_size/15),  
            x = xloc, 
            y = yloc,  
            yanchor="top", 
            title = var
        )
        zmin = 0
        zmax = len(bvals)-1
    else:
        cbar = dict(
            len=continuous_height,
            x = xloc,
            y=yloc,
            yanchor="top",
            title = var
        )
        zmin = None
        zmax = None
    
    hover_text = xcol +': %{x}<br>'+                 _get_custom(meta[[var]])+                 '<br>'
    
    hm = go.Heatmap(
        name = var, 
        z = meta[varcol].values, 
        x = meta[xcol].values, 
        y = yvals,
        zmin = zmin,
        zmax = zmax,
        customdata = meta[[var]],
        # hovertext=meta[var].values,
        hovertemplate= hover_text,
        colorscale = colorscale,
        colorbar = cbar,
        legendgroup = legendgrp,
        legendgrouptitle = {'text':var}
    )
    
    return(hm)


def trim_quantile(x, q, direction='gt'):
    qval = np.quantile(x, q) 
    res = x.copy()
    if direction == 'gt':
        res[res>qval] = qval
    if direction == 'lt':
        res[res<qval] = qval
    return(res)


### Annotated DEG Heatmap

def python_match(orderedvals, values):
    '''Get indexes of values matching to a set of ordered values
    
    Parameters
    ----------
    orderedvals : list
        list of ordered values, expected to be unique.
    values : list
        list of all values to order
    
    Returns
    -------
    list
        Integer list of indexes in values, ordered by orderedvals.
    
    Examples
    --------
    
    
    '''
    imatch = [[i for i in range(len(values)) if values[i] == y] for y in orderedvals]
    imatch = [item for sublist in imatch for item in sublist] # unnest

    return(imatch)


# %%


def map_p_binary(pvals, alpha=0.05, type_dict = {'s':1,'ns':-1, 'other':0}):
    ''' Map p-values to binary symbols or values based on a single significance thresholds
    '''
    
    res = [type_dict['s'] if x <= alpha else type_dict['ns'] if x > alpha else type_dict['other'] for x in pvals]

    return(res)

# mynum = [0,0.049,0.05,0.051,1,float("NaN")]
# pdict_num = {'s':1,'ns':-1, 'other':0}
# map_p_binary(mynum, type_dict = pdict_num)

# pdict_sym = {'s':"*",'ns':"", 'other':""}
# map_p_binary(mynum, type_dict = pdict_sym)


def differential_hm(
    df, 
    y_col = 'gene',
    obs_col = 'obs', 
    # feat_col = 'gene', 
    signif_col = 'adjP', 
    es_col='logFC', 
    filter_dict = dict(),
    signif_cutoff = 0.05,
    cluster_rows = False, 
    cluster_cols = False,
    reverse_rows = False, 
    meta_plot = [], 
    meta_plot_order = None,
    symbol_size = 10,
    verbose = False, 
    meta_colors = 'viridis',
    color_scale_hm= 'RdBu_r',
    title = "DEG Heatmap",
    # drop_missing=True,
    fontsize = 12, 
    legend_font_size = 10, 
    legend_column_spacing = 0.1, 
    legend_row_spacing = 0.01,
    legend_height_per_level = 0.03, 
    legend_continuous_height = 0.2, 
    legend_matrix_val_height = 0.2,
    check = False
    ):
    df_raw = df.copy()

    # Filter data (metadata or features okay)
    df_fmt = df_raw.copy()
    for key in filter_dict:
        df_fmt = df_fmt.loc[df_fmt[key].isin(filter_dict[key])]
    
    # metadata
    if meta_plot is not None:
        # format colors
        if not isinstance(meta_colors, dict):
            meta_colors = pd.Series(meta_colors, index = meta_plot).to_dict()
        else:
            missing_cols=[x for x in meta_plot if x not in meta_colors.keys()]
            for mvar in missing_cols:
                meta_colors.update({mvar : 'viridis'})
        keep_cols = meta_plot.copy()
        if obs_col not in keep_cols:
            keep_cols.append(obs_col)
    else:
        keep_cols = obs_col
    meta = df_fmt[keep_cols].copy()
    meta = meta.drop_duplicates()
    if check:
        print('Metadata raw')   
        print(meta[obs_col].values)   
    
    # format meta after filter
    meta = reset_all_categorical(meta)
    if meta_plot_order is not None:
        meta = meta.sort_values(by = meta_plot_order, na_position = 'last')
        meta.reset_index(drop =True)
    if check:
        print('Metadata sorted')   
        print(meta[obs_col].values)   
        print('column order: ')
        print(meta_plot_order)
            
    # convert es values to matrix
    mat_es = df_fmt.copy().pivot(index=y_col, columns=obs_col, values = es_col)
    
    # pvalue recoded matrices
    pdict_sym = {'s':"*",'ns':"", 'other':""}
    pdict_num = {'s':1,'ns':-1, 'other':0}
    
    mat_signif = df_fmt.copy().pivot(index=y_col, columns=obs_col, values = signif_col)
    
    mat_signif_text = mat_signif.copy().apply(lambda x: map_p_binary(x, alpha = signif_cutoff, type_dict = pdict_sym), axis = 1, result_type='expand')
    mat_signif_text = mat_signif_text.rename(columns={i:mat_signif.columns[i] for i in range(mat_signif.shape[1])}, errors="raise")

    # order columns of matrixes
    if cluster_cols:
        if verbose:
            print("clustering columns")
            print("Matrix shape: {}".format(mat_es.shape))
        col_clorder_index = cluster_matrix(mat_es, index = 1, fill_na = 0, value_return='index')
        imeta = python_match(orderedvals=col_clorder_index, values=meta[obs_col].values.tolist())
        meta = meta.iloc[imeta,]
        meta = meta.reset_index(drop=True)
        mat_signif_text = mat_signif_text.loc[:, col_clorder_index]
        mat_signif = mat_signif.loc[:, col_clorder_index]
        mat_es = mat_es.loc[:, col_clorder_index]
    else:
        mat_signif = mat_signif.loc[:, meta[obs_col].values.tolist()]
        mat_signif_text = mat_signif_text.loc[:, meta[obs_col].values.tolist()]  # need tolist() to get strings instead of categories.
        mat_es = mat_es.loc[:, meta[obs_col].values.tolist()]
    if check:
        print("pval matrix after column ordering")
        print(mat_signif)
        print(meta[obs_col].values)
    
    # order rows
    if cluster_rows:    
        if verbose:
            print("clustering rows")
            print("Matrix shape: {}".format(mat.shape))
        row_clorder_index = cluster_matrix(mat_es, index = 0, fill_na = 0, value_return='index')
        mat_signif = mat_signif.loc[row_clorder_index, :]
        mat_signif_text = mat_signif_text.loc[row_clorder_index, :]
        mat_es = mat_es.loc[row_clorder_index, :]
    elif reverse_rows:
        if verbose:
            print("reversing rows")
            print("Matrix shape: {}".format(mat.shape))
        mat_signif = mat_signif.loc[::-1]
        mat_signif_text = mat_signif_text.loc[::-1]
        mat_es = mat_es.loc[::-1]
        
    if check:
        print("pval matrix after row ordering")
        print(mat_signif)
        print(meta[obs_col].values)
    
    # Construct Plot
    
    # hover pattern
    hover_text_mat = obs_col +': %{x}<br>'+                     y_col +': %{y}<br>'+                     es_col +': %{z:.3f}<br>'+                     signif_col+': %{customdata:.3f}<br>'+                     'significance: %{text}<br>'
    
    # plot
    if meta_plot is not None: # metadata annotations
        # initialize a subplot matrix
        nrows =  len(df_fmt[y_col].unique())
        rheights = _calc_subplot_row_heights(nfeat = nrows, nmeta = len(meta_plot))
        fig = plotly.subplots.make_subplots(
            rows=1+len(meta_plot), 
            cols=1, 
            row_heights = rheights, 
            vertical_spacing = 0.01,
            shared_xaxes=True
        )   

        # Calculate subplot legend positions
        xloc, yloc = _calc_meta_legend_locs(
            meta = meta, 
            meta_plot = meta_plot, 
            height_per_level = legend_height_per_level, 
            continuous_height = legend_continuous_height, 
            matrix_val_height = 0.2, 
            legend_font_size = legend_font_size, 
            xshift = legend_column_spacing, 
            yshift = legend_row_spacing
        )

        # Metadata annotation plots
        for i in range(len(meta_plot)):
            var = meta_plot[i]
            if verbose:
                print('making hm trace for metadata {}'.format(var))
            hm_trace = format_heatmap_anno(
                meta = meta, 
                var=var, 
                xcol=obs_col,
                colorscale=meta_colors[var],
                legendgrp= str(i), 
                xloc = xloc[i], 
                yloc = yloc[i], 
                legend_font_size=legend_font_size)

            fig.append_trace(
                hm_trace,
                row=i+1, 
                col=1
            )
            
        # Matrix plot
        cbar_mat = dict(
            len=legend_matrix_val_height,
            y=yloc[len(meta_plot)],
            yanchor='top',
            x = xloc[len(meta_plot)],
            title = title,
        )
        dat_hm = go.Heatmap(
            z = mat_es, 
            y = mat_es.index, 
            x = mat_es.columns, 
            customdata = mat_signif,
            hovertemplate=hover_text_mat,
            text = mat_signif_text,
            texttemplate="%{text}",
            textfont={"size":20},
            colorbar = cbar_mat, 
            colorscale=color_scale_hm,
            zmid=0,
            legendgroup = str(len(meta_plot)+1),
            name = title)

        fig.append_trace(
            dat_hm,
            row=len(meta_plot)+1, 
            col=1
        )        
    else:
        
        fig = go.Heatmap(
            z = mat_es, 
            y = mat_es.index, 
            x = mat_es.columns, 
            customdata = mat_signif,
            hovertemplate=hover_text_mat,
            text = mat_signif_text,
            texttemplate="%{text}",
            textfont={"size":symbol_size},
            colorbar = dict(
                len=legend_matrix_val_height,
                y=1,
                yanchor='top',
                x = 1,
                title = title,
            ), 
            legendgroup = 1,
            name = title)
    
    return(fig)


def cluster_matrix(mat, index=0, fill_na = 0, value_return = 'order'):
    '''
    Cluster a matrix by row (index = 0) or column (index = 1) using hierarchical
    clustering (sklearn AgglomerativeClustering). Return the sorted column order
    as numeric locations (value_return = 'order'
    
    Examples
    --------
    ncells = 100
    nfeat = 10
    meta_levels = {
        'group' : ['A','B','C'],
        'timepoint' : [1,2,3,4],
        'fruit' : ['apple','banana']
        }
    meta = pd.DataFrame([])
    for key in meta_levels:
        meta[key] = pd.Series([random.sample(meta_levels[key], 1)[0] for i in range(ncells)], dtype="category")
    meta['cell'] = ['cell_' + str(x) for x in range(ncells)]
    mat = np.array([np.random.negative_binomial(2, 0.8, size=ncells)  for i in range(nfeat)])

    feat = ["feat"+str(i+1) for i in range(nfeat)]

    matdf = pd.DataFrame(mat)
    matdf['feat'] = feat
    matdf = matdf.set_index('feat')
    
    irow = cluster_matrix(matdf)
    matdf.iloc[irow,]
    
    icol = cluster_matrix(matdf, index = 1)
    matdf.iloc[:,icol]
    
    '''
    if index == 0:
        mat_cl = mat.copy()
    elif index == 1:
        mat_cl = mat.transpose().copy()
    else:
        print("index should be 0 for clustering rows or 1 for clustering columns")
        
    mat_cl = mat_cl.fillna(fill_na)
    
    cl_res_list = []

    for i in range(mat_cl.shape[0]):
        hclust_col = AgglomerativeClustering(n_clusters = i+1, metric = 'euclidean', linkage='ward')
        cname = 'itr'+str(i+1)
        cl_res_list.append(pd.DataFrame({cname:hclust_col.fit_predict(mat_cl)}))
        
    cl_res = pd.concat(cl_res_list, axis=1)
        
    dim_order = cl_res.sort_values(by = cl_res.columns.tolist()).index
    
    if value_return == 'order':
        return(dim_order.tolist())
    elif value_return == 'index':
        return(mat_cl.index[dim_order].tolist())


### Parsing

# used to parse pasted in input genes
def parse_input(in_string, delim = [',','"',"'"]):
    for sym in delim:
        in_string = re.sub(sym,"\t", in_string)
    in_string = re.sub("\\s+","\t", in_string)
    out_string = in_string.split('\t')
    out_string = [x for x in out_string if x != '']
    return(out_string)

#' Return input values case-matched to a reference 'valid list'
#' for all value matches found, otherwise otherwise return original value
def match_valid_nocase(in_list, valid_list):
    in_list_lower = [x.lower() for x in in_list]
    valid_list_lower = [x.lower() for x in valid_list]
    imatch = [i for i in range(len(in_list)) if in_list_lower[i] in valid_list_lower]
    out_list = in_list.copy()
    for x in imatch:
        out_list[x] = valid_list[valid_list_lower.index(in_list_lower[x])]
    return out_list


### Data formatting
def check_cat_dict(cat_dict, df, verbose = False):
    ''' 
    Checks that a dictionary of column names and ordered expected values represents complete 
    data given the data frame.
    Requires: from warnings import warn
    '''
    
    for name in cat_dict.keys():
        if(name not in df.columns):
            warn('category dictionary column name not found in data: ' + name)
        else:
            data_cat_all = pd.unique(df[name])
            missing_cat = [x for x in data_cat_all if x not in cat_dict[name]]
            if len(missing_cat) > 0:
                warn(
                    'The following values in data column "' + name + '" are not found in category dictionary: ' + ", ".join(missing_cat) + 
                    ". Add to dictionary if you wish to include as categorical levels for plot filtering, etc."
                )
            extra_cat = [x for x in cat_dict[name] if x not in data_cat_all]
            if len(extra_cat) > 0:
                warn(
                    'Extra data levels found in data dictionary for column "' + name + '": ' + ", ".join(extra_cat) + 
                    ". These will be removed as categories during variable formatting."
                )
    if(verbose):
        print("Finished checking category dictionary")

def format_df(df, cat_dict={}):
    res = df.copy()
    
    # Formats
    types = [res[x].dtype for x in res.columns]
    
    # Convert string to categorical
    i_str = [i for i in range(len(types)) if types[i] =='O']
    if(len(i_str) > 0):
        if(len(i_str) > 0):
            cnames = [res.columns[i] for i in i_str]
            cat_dict_all = {x:cat_dict[x] if x in cat_dict.keys() else np.sort(pd.unique(res[x])) for x in cnames}
            for col in cnames:
                res[col] = res[col].astype(pd.CategoricalDtype(categories=cat_dict_all[col], ordered = True))
    return res

def reset_filtered_categorical(cat_val):
    '''reset a categorical datatype to only the values represented in the dataset
    
    Params
    ------
    cat_val Categorical data series
        For example, extracted from a pandas dataframe df[cat_val]
    
    '''
    orig_cat = cat_val.dtype.categories.values
    subset_cat = [x for x in orig_cat if x in cat_val.values]
    subset_dtype = pd.CategoricalDtype(subset_cat, ordered=True)
    new_val = cat_val.astype(subset_dtype)
    
    return(new_val)

def reset_all_categorical(df):
    '''reset all categorical datatype to only the values represented in the dataset
    
    Params
    ------
    cat_val Categorical data series
        For example, extracted from a pandas dataframe df[cat_val]
    
    '''
    df_relevel = df.copy()
    for var in df_relevel.columns:
        if isinstance(df_relevel[var].dtype, pd.CategoricalDtype):
            df_relevel[var] = reset_filtered_categorical(df_relevel[var])
    
    return(df_relevel)

## Data Prep-- Modify for new datasets

### Constants
app_title = 'VRd In-Vitro T-cell DEGs'

no_selection_instruction = 'Please make a gene selection to view gene information'

hm_defaults = ['comparison','treatment']
volcano_ct_default = 'CD4 Naive'
# volcano_experiment_default = 'All Subjects'
volcano_comparison_default = 'TEA-seq Dex. 4 hr'

debug = False

### DEG Results
fp = './data/all_mast_deg_2023-09-06.parquet'
df = pl.read_parquet(fp).to_pandas()

feat = df.gene.unique().tolist()

#### Format Categorical Data

ct_order = ['CD4 Naive','CD4 CM','CD4 EM','CD4 Treg',
            'CD8 Naive', 'CD8 Memory']
expt_order = ['VRd Perturbations']
comp_order = [
    'TEA-seq Bor. 4 hr','TEA-seq Bor. 24 hr','TEA-seq Bor. 72 hr',
    'TEA-seq Len. 4 hr','TEA-seq Len. 24 hr','TEA-seq Len. 72 hr',
    'TEA-seq Dex. 4 hr','TEA-seq Dex. 24 hr'
]
treat_order = ['VRd','bortezomib','lenalidomide','dexamethasone']

# use this to populate possible volcano plot options
deg_result_combos = df.loc[:,['cell_type','comparison']].drop_duplicates()

# Make Categorical data dictionary to explicitly level categorical data. run the `check_cat_dict()` function to check for potential errors


#### Colors
fp_col = './data/meta_color_dict.json'
with open(fp_col, 'r') as handle:
    preset_colors_dict = json.load(handle)

### Pathway genes
fp_gs = './data/custom_gs_dict.json'
with open(fp_gs, 'r') as handle:
    all_gs = json.load(handle)

## App
pio.templates.default = 'plotly_white'

# ### Selector Panel

selector_panel_row = dbc.Row(
    className = 'p-4',
    id='selector_panel',
    children = [
        dbc.Col(
            width = 2, 
            children = [
                html.Div(
                    id='volcano_comp_dropdown_container',
                    children=[
                        html.P('DEG Comparison'),
                        dcc.Dropdown(
                            options = [{'label': i, 'value': i } for i in comp_order], 
                            id ='volcano_comp_dropdown', 
                            multi= False, 
                            clearable=False,
                            value = volcano_comparison_default
                        )
                    ]
                ),
            ]
        ),
        dbc.Col(
            width = 2, 
            children = [
                html.Div(
                    id='volcano_ct_dropdown_container',
                    children=[     
                        html.P('Cell Type'),
                        dcc.Dropdown(
                            id ='volcano_ct_dropdown', 
                            multi= False, 
                            clearable=False,
                        )
                    ]
                )
            ]
        ),
        dbc.Col(
            width = 2, 
            children = [
                html.Div(
                    id='p_input_container',
                    children=[
                        html.P('Select Significance Cutoff'),
                        dcc.Dropdown(
                            options = [{'label': i, 'value': i } for i in [0.01,0.05,0.1]], 
                            id ='p_input', 
                            multi= False, 
                            clearable=False,
                            value = 0.05
                        )
                    ]
                )
            ]
        ),
        dbc.Col(
            width = 2, 
            children = [
                html.Div(
                    id='es_input_container',
                    children=[     
                        html.P('Select Effect Size Cutoff'),
                        dcc.Dropdown(
                            options = [{'label': i, 'value': i } for i in [0, 0.1, 0.2, 0.5, 1]], 
                            id ='es_input', 
                            multi= False, 
                            clearable=False,
                            value = 0.1
                        )
                    ]
                )
            ]
        ),
        dbc.Col(
            width = 2, 
            children = [
                html.Div(
                    id='gene_hm_container', 
                    children=[
                        html.P('Heatmap Gene Selection'),
                        dcc.Dropdown(
                            options = [{'label': i, 'value': i } for i in feat], 
                            id ='select_hm_gene', 
                            multi= False, 
                            value = None,
                            placeholder="Click from volcano or select a gene")
                    ]
                ),  
            ]
        )
    ]
)

## Volcano details 
volcano_annotation_div = dbc.Card(
    id='master_gene_selection_div',
    children = [
        dbc.CardHeader([
            html.H2(
                "Select Genes for Volcano Plot Annotation",
                style = {
                    'color':colors["--bs-hise-blue-1"],
                    'fontSize': '20px'
                }
            ),
            html.P(
                "Use any combination of input methods to select genes to annotate",
                style = {'color':colors["--bs-hise-blue-1"]}
            )
        ]),
        dbc.CardBody([
            dbc.Row(dbc.Col(
                html.Div(
                    id='gene_selector_container',
                    children=[
                        html.B('Select from Dropdown:'),
                        html.Div(id='dropdown_col', 
                            children=[     
                                dcc.Dropdown(
                                    options = [{'label': i, 'value': i } for i in feat], 
                                    id ='select_gene', 
                                    multi= True, 
                                    value = None,
                                    clearable = True,
                                    placeholder="Select gene(s) from dropdown")
                            ]),
                    ]
                )
            )),
            dbc.Row(dbc.Col(
                html.Div(
                    id='manual_input_container',
                    children=[
                        html.B('Paste Gene List:'),
                        html.Div(
                            id='manual_input_column',
                            children=[     
                                dcc.Textarea(
                                    id ='manual_input', 
                                    value = None,
                                    style = {'width': '100%'},
                                    placeholder="Paste list of gene names (comma-, white-space-, or quotation-delimited)"
                                ),
                                html.P(
                                    id = 'manual_input_notes',
                                    style = {
                                        'size': 10, 
                                        'color':"gray",
                                        'overflow-y': 'auto', 
                                        'maxHeight':'50px'
                                    },
                                )
                            ]
                        )
                    ]
                )
            )),
            dbc.Row(dbc.Col(
                html.Div(
                    id='pathway_selector_container',
                    children=[
                        html.B('Select Gene Set:'),
                        html.Div(
                            id='dropdown_col_pw', 
                            children=[     
                                dcc.Dropdown(
                                    options = [{'label': key, 'value': key, 'title': all_gs[key]['description']} for key in all_gs], 
                                    id ='select_pathway', 
                                    multi= False, 
                                    value = None,
                                    placeholder="Select a geneset"
                                ),
                                html.P(
                                    id = 'select_pathway_notes',
                                    style = {
                                        'size': 10, 
                                        'color':"gray",
                                        'overflow-y': 'auto', 
                                        'maxHeight':'50px'
                                    },              
                                )
                            ]
                        )
                    ]
                )
            )),
            dbc.Row([
                html.Div(
                    id='clicked_data_container',
                    children=[
                        html.B('Plot-selected Genes:'),
                        html.Div(
                            id='dropdown_col_clickgenes', 
                            children=[
                                dbc.Row([
                                    dbc.Col(
                                        width = 10,
                                        children = dcc.Dropdown(
                                            options = [{'label': i, 'value': i } for i in feat], 
                                            id ='dropdown_plot_selected', 
                                            multi= True, 
                                            value = None,
                                            clearable = False,
                                            disabled = True,
                                            style = {'maxHeight': '150px', 'overflow-y': 'auto'},
                                            placeholder="Selected gene(s) from volcano plot",
                                            className = 'd-inline'
                                        )
                                    ),
                                    dbc.Col(
                                        width = 2,
                                        children = dbc.Button(
                                            n_clicks = 0,
                                            id = 'clear_plot_selection_button',
                                            title = "Clear all markers selected from plot",
                                            color = 'danger',
                                            children = [
                                                DashIconify(
                                                    icon = 'fa:close',
                                                    height = 20
                                                )
                                            ],
                                            size = 'sm',
                                            className = 'd-inline'
                                        )
                                    )
                                ])
                            ]
                        )
                    ]
                )
            ]),
            html.Hr(),
            dbc.Row(dbc.Col(
                children = [
                    html.B(
                        'Selections:'
                    ),
                    html.Div(
                        id='view_gene_sel_container',
                        children=[
                            dcc.Dropdown(
                                id='select_gene_view',
                                className = 'col-12',
                                options = [{'label': i, 'value': i } for i in feat],
                                clearable = True,
                                multi = True,
                                value = None,
                                style = {'maxHeight': '100px', 'overflow-y': 'auto'},
                                placeholder = 'Make selections above to build a gene annotation list',
                                disabled = True
                            ),
                        ]
                    )
                ]),
            ),
            html.Br(),
            dbc.Row(dbc.Col(
                html.Div(
                    id='gene_button_row', 
                    className = 'text-center',
                    children=[
                        dbc.Button(
                            id = 'update-genes',
                            title = "Update Gene Annotation",
                            className = 'bodyButton col-11',
                            n_clicks=0,
                            outline = False,
                            children = [
                                'Update Plot',
                                DashIconify(
                                    id = 'volcano_gene_button_icon',
                                    icon = 'fa:arrows-rotate',
                                    height = 20,
                                    color = 'black'
                                )                  
                            ]
                        )
                    ]
                )
            ))
        ])
    ]
)

## Heatmap details 
hm_details_container_div = html.Div(
    id = 'hm_details_container',
    children = [
        html.H2(
            "Customize Heatmap",
            style = {
                'color':colors["--bs-hise-blue-1"],
                'fontSize': '20px'
            }
        ),
        dbc.Row(
            id='hm_select_r1',
            children = [
                dbc.Col(
                    id='metadata_plot_dropdown', 
                    width = '3',
                    children=[
                        html.P('Select Metadata Column(s) to Plot', className='instruction1'),
                        dcc.Dropdown(
                            # options = [{'label': i, 'value': i } for i in df.select_dtypes(include=['category']).columns.values], 
                            options = hm_defaults,
                            id ='select_hm_meta', 
                            multi= True, 
                            value = hm_defaults
                        )
                    ]
                ),
                dbc.Col(
                    id='metadata_sort_checklist', 
                    width = 2,
                    children=[     
                        html.P('Sort cells by metadata?', className='instruction1'),
                        dcc.Checklist(
                            id ='select_hm_sort',
                            inputStyle={"margin-left": "10px","margin-right": "2px"}, 
                            options=[],
                            value=[]
                        )
                    ]
                ),
                dbc.Col(
                    id='cluster_col_radio_container', 
                    width = 2,
                    children=[     
                        html.P('Cluster columns?', className='instruction1'),
                        dcc.RadioItems(
                            id ='cluster_col_radio', 
                            options=['No', 'Yes'],
                            value='No',
                            inputStyle={"margin-left": "10px","margin-right": "2px"}
                        )
                    ]
                ),
                dbc.Col(
                    id='cluster_row_radio_container', 
                    width = 2,
                    children=[     
                        html.P('Cluster rows?', className='instruction1'),
                        dcc.RadioItems(
                            id ='cluster_row_radio', 
                            options=['No', 'Yes'],
                            value='No',
                            inputStyle={"margin-left": "10px","margin-right": "2px"}
                        )
                    ]
                ),
            ]
        ),
    ]
)

# ### Plot Panel

plot_row = dbc.Row(
    className = 'px-4',
    id='plot_panel',
    children=[
        dbc.Col(
            id='volcano_column', 
            width = 5,
            style = {
                'borderRight': '2px solid #003056'
            },
            children = [
                dbc.Row(
                    className = 'px-1',
                    children = [
                        dbc.Col(
                            dbc.Button(
                                id = 'volcano_annotation_button',
                                n_clicks=0,
                                outline = False,
                                children = [
                                    DashIconify(
                                        id = 'volcano_annotation_button_icon1',
                                        icon = "fa:pencil"
                                    ), ' ',
                                    DashIconify(
                                        id = 'volcano_annotation_button_icon2',
                                        icon = "fa:caret-down"
                                    )
                                ]
                            )
                        )
                    ]
                ),
                dbc.Row(className = 'p-2', children = volcano_annotation_div),
                dbc.Row([
                    html.Div(
                        id='volcano_plot-Container',
                        children = dcc.Graph(
                            id = 'volcano_plot', 
                            responsive = True
                        )
                    )
                ])
            ]
        ),
        dbc.Col(
            id='hm_column', 
            width = 7,
            children = [
                dbc.Row(
                    className = 'px-1',
                    children = [
                        dbc.Col(
                            dbc.Button(
                                id = 'hm_details_button',
                                n_clicks=0,
                                outline = False,
                                children = [
                                    DashIconify(
                                        id = 'hm_details_button_icon1',
                                        icon = 'fa:gear',
                                        height = 20
                                    ),
                                    DashIconify(
                                        id = 'hm_details_button_icon2',
                                        icon = 'fa:caret-down',
                                        height = 20
                                    )
                                ]
                            )
                        )
                    ]
                ),
                dbc.Row(
                    className = 'p-2',
                    children = hm_details_container_div
                ),
                dbc.Row(
                    html.Div(id = 'heatmap-container')
                )
            ]
        )
    ]
)

### Gene Panel

gene_row = dbc.Row(
    id = 'gene-info-row',
    className = 'p-4',
    children = [
        html.Div(
            id = 'gene-info-container',
            children = []
        )
    ] 
)

# ### Layout

app.config.suppress_callback_exceptions=True

app_content = html.Div(
    children=[
    
        dcc.Store(id='gene_select_store', storage_type='memory'),
        dcc.Store(id='volcano_click_select_store', storage_type='memory'),
        dcc.Store(id='hm_gene_store', storage_type='memory'),

        selector_panel_row,
        plot_row,
        gene_row,
    
])# end layout

# ### Callback-app

# CallBacks ===================================================================================#
# Header

@app.callback(
    Output('popover', 'is_open'),
    [Input('info_button', 'n_clicks')],
    [State('popover', 'is_open')],
)
def update_callback(n, is_open):
    if n :
        return not is_open
    else :
        return is_open

@app.callback(
    Output('manual_input_notes', 'children'),
    Input('manual_input', 'value'))
def update_man_input(man_vals):
    if man_vals is None or man_vals == "":
        outstr = ""
    else:
        man_vals = parse_input(man_vals)
        man_vals = match_valid_nocase(man_vals, feat)  # match case to expected if found in feature list
        missing_vals = [x for x in man_vals if x not in feat]
        if len(missing_vals) > 0:
            missing_str = ", ".join(missing_vals)
            outstr = 'The following input values are not present in dataset: "{}"'.format(missing_str)
        else:
            outstr = 'All values detected in dataset'
    return outstr

@app.callback(
    Output('select_pathway_notes', 'children'),
    Input('select_pathway', 'value'))
def update_pw_notes(pw):
    if pw is None:
        outstr = ""
    else:
        pwvals = all_gs[pw]['features']
        missing_vals = [x for x in pwvals if x not in feat]
        if len(missing_vals) > 0:
            missing_str = ", ".join(missing_vals)
            outstr = 'The following geneset features are not present in dataset: "{}"'.format(missing_str)
        else:
            outstr = 'All values detected in dataset'
    return outstr

@app.callback(
    Output('gene_select_store', 'data'),
    [Input('update-genes', 'n_clicks')],
    [State('select_gene_view', 'value')])
def update_gene_store(n, genes):
    if n is not None:
        return genes
    else:
        return None
    
# button click
@app.callback(
    [Output('master_gene_selection_div', 'className'),
     Output('volcano_annotation_button', 'title'),
     Output('volcano_annotation_button_icon1', 'icon'),
     Output('volcano_annotation_button_icon2', 'icon')],
    [Input('volcano_annotation_button', 'n_clicks')]
)
def button_toggle(n_clicks):
    if (n_clicks is None) or (n_clicks%2 == 0) :
        return ['noshow', "Add Volcano Plot Annotations", 'fa:pencil','fa:caret-down']
    else :
        return ['expandPanel', 'Close Volcano Plot Annotation Panel', 'fa:pencil','fa:caret-up']

# ### Callback-volcano

# Tab 1, Volcano Panel
@app.callback(
    Output('volcano_ct_dropdown', 'options'),
    Output('volcano_ct_dropdown', 'value'),
    Input('volcano_comp_dropdown','value'))
def update_volcano_ct(comp):
    if comp is not None and len(comp) >0:
        ct_options = deg_result_combos.sort_values('cell_type').loc[deg_result_combos['comparison'] == comp,'cell_type'].drop_duplicates().tolist()
        outopt = [{'label': x, 'value': x} for x in ct_options]
        outval = ct_options[0]
    else:
        outopt = []
        outval = None
    return outopt, outval

@app.callback(
    [Output('volcano_plot', 'figure'),
     Output('volcano_plot', 'config')],
    [Input('gene_select_store', 'data'),
    Input('p_input', 'value'),
    Input('es_input', 'value'),
    Input('volcano_ct_dropdown', 'value'),
    Input('volcano_comp_dropdown', 'value')]
)
def update_volcano(gene, p, es, ct, comp):
    if debug == True:
        print("celltype: "+ ct)
        print("comparison: "+ comp)
    # Make Plot list from input deg df's
    # fig_ls = []
    temp = df.loc[(df.cell_type == ct) & (df.comparison == comp),:].copy()
    
    # # add change variable
    # temp
    
    plt = build_volcano(
        p_name='adjP', 
        es_name='logFC',
        hovername = 'gene',
        highlight = gene,
        de_df = temp, 
        es_cutoff = es, 
        alpha = p,
        opacity_val=0.5,
        title = comp + "<br>"+ct,
        title_size = 16,
        height = 500
    )
    plt.update_layout(
        autosize = True,
        margin=dict(
                l=20,
                r=20,
                b=25,
                t=75,
                pad=5
        )
    )
    
    cf = {
        'toImageButtonOptions':{
            'filename': 'degapp_volcano_'+str(comp)+'_'+str(ct)+'_alpha' + str(p)+'_es'+str(es),
            'format': 'svg'
        }
    }
    
    return [plt, cf]

@app.callback(
    Output('dropdown_plot_selected', 'value'),
    Input('volcano_click_select_store', 'data'))
def update_select_view(click_list):
    return click_list

@app.callback(
    Output('select_gene_view', 'value'),
    [Input('volcano_click_select_store', 'data'),  # list of clickData inputs for all volcano plots
    Input('select_pathway', 'value'),
    Input('manual_input', 'value'),
    Input('select_gene', 'value')])
def update_gene_dropdown(click_list, pw, man_vals, currentvals):
    if debug == True:
        print("pw: ")
        print(pw)
        print("man_vals:")
        print(man_vals)
        print("clicklist:")
        print(click_list)
        print('currentvalues:')
        print(currentvals)
    if currentvals is None:
        currentvals = []
    if pw is not None:
        pwvals = all_gs[pw]['features']
        for pwval in pwvals:
            if (pwval not in currentvals) and (pwval in feat):
                currentvals.append(pwval)
    if man_vals is not None:
        man_vals = parse_input(man_vals)
        man_vals = match_valid_nocase(man_vals, feat)  # match case to expected if found in feature list
        for manval in man_vals:
            if (manval not in currentvals) and (manval in feat):
                currentvals.append(manval)
    for clickData in click_list:
        if clickData is not None:
            if (clickData not in currentvals) and (clickData in feat):
                currentvals.append(clickData)
    return currentvals

@app.callback(
    [Output('volcano_click_select_store', 'data'),
    Output('clear_plot_selection_button', 'n_clicks')],
    [Input('clear_plot_selection_button', 'n_clicks'),
     [Input('volcano_plot', 'clickData')],
     [Input('volcano_plot', 'selectedData')],
     State('volcano_click_select_store', 'data')]
  )
def update_gene_dropdown(n_clicks, click_list, select_list, currentvals):
    if(n_clicks > 0):
        currentvals = []
    else:
        if currentvals is None:
            currentvals = []
        for clickData in click_list:
            if clickData is not None:
                addval = clickData['points'][0]['hovertext']
                if (addval not in currentvals):
                    currentvals.append(addval)
        for selectData in select_list:
            if selectData is not None:
                for i in range(len(selectData['points'])):
                    addval = selectData['points'][i]['hovertext']
                    if (addval not in currentvals) and (addval in feat):
                        currentvals.append(addval)
    return [currentvals, 0]

# ### Callback-HM

# Tab 3
@app.callback(
    Output('hm_gene_store', 'data'),
    Input('select_hm_gene', 'value'))
def update_hmgene_store(gene):
    if gene:
        return gene

@app.callback(
    Output('select_hm_gene', 'value'),
    [Input('volcano_plot', 'clickData')])  
def update_gene_dropdown(click_list):
    if debug == True:
        print(click_list)

    if click_list is None:
        return 'CXCR4'
    
    val = click_list['points'][0]['hovertext']
    
    return val

@app.callback(
    [Output('select_hm_sort', 'options'),
    Output('select_hm_sort', 'value')],
    [Input('select_hm_meta','value')])
def update_hm_metadata(meta_val):
    if meta_val is not None and len(meta_val) > 0:
        if type(meta_val) != list:
            meta_val = [meta_val]
        outval = [{'label': x, 'value': x} for x in meta_val]
        selection = [x for x in meta_val]
        
    else:
        outval = []
        selection = []
    return outval, selection
    
@app.callback(
    [Output('select_hm_filter_value', 'options'),
    Output('select_hm_filter_value', 'value')],
    [Input('select_hm_filter_col','value'),
    Input('volcano_ct_dropdown', 'value'),
    Input('volcano_comp_dropdown', 'value')])
def update_hm_filter(meta_val, ct, comp):
    if meta_val is not None:
        outopt = [{'label': x, 'value': x} for x in df[meta_val].unique()]
        if meta_val == 'cell_type':
            outval = ct
        elif meta_val == 'comparison':
            outval = comp
    else:
        outval = None
        outopt = []
    return outopt, outval

@app.callback(
    [Output('hm_details_container', 'className'),
     Output('hm_details_button', 'title'),
     Output('hm_details_button_icon1', 'icon'),
     Output('hm_details_button_icon2', 'icon')],
    [Input('hm_details_button', 'n_clicks')]
)
def button_toggle(n_clicks):
    if (n_clicks is None) or (n_clicks%2 == 0) :
        return ['noshow', "Heatmap Options", 'fa:gear','fa:caret-down']
    else :
        return ['expandPanel', 'Close Heatmap Options', 'fa:gear','fa:caret-up']

    
@app.callback(
    Output('heatmap-container', 'children'),
    [Input('hm_gene_store', 'data'),
     Input('p_input', 'value'),
     Input('select_hm_meta', 'value'),
     Input('select_hm_sort', 'value'),
     Input('cluster_col_radio', 'value'),
     Input('cluster_row_radio', 'value')]
)
def update_heatmap(value, signif_cutoff, metaplot, metasort, cluster_col, cluster_rows):
    if value is None:
        outfig = [html.Div("Make a gene selection to view DEG heatmap for all comparisons")]
    elif (type(value) == list) & (value[0] is None):
        outfig = [html.Div("Make a gene selection to view DEG heatmap for all comparisons")]
    else:
        if metasort is None or len(metasort)<1 or metaplot is None or len(metaplot) <1:
            meta_sort = None
        else:
            meta_sort = [x for x in metasort if x in metaplot]

        # translate clustering input
        cluster_dict={'Yes':True, 'No':False}
        cluster_col_b = cluster_dict[cluster_col]
        cluster_rows_b = cluster_dict[cluster_rows]

        if metaplot is None or len(metaplot) <1:
            color_dict_plot = 'plasma'
        else:
            color_dict_plot = dict()
            for val in metaplot:
                if val in preset_colors_dict.keys():
                    color_dict_plot.update({val:preset_colors_dict[val] })
                else:
                    color_dict_plot.update({val:'plasma'})
                    
        # reverse the celltype order for row plot

        if type(value)==str:
            value = [value]            
        fdict={'gene':value}

        fig = differential_hm(
            df=df, 
            filter_dict=fdict,
            y_col = 'cell_type',
            reverse_rows = True,
            cluster_cols = cluster_col_b, 
            cluster_rows = cluster_rows_b,
            signif_cutoff = signif_cutoff,
            signif_col = 'adjP', 
            es_col = 'logFC', 
            obs_col = 'comparison', 
            meta_plot = metaplot,
            meta_plot_order=meta_sort,
            meta_colors=color_dict_plot,
            fontsize = 12,
            legend_font_size = 10,
            legend_row_spacing = 0.005,
            legend_column_spacing=0.4
        )
        fig.update_layout(
            autosize = True,
            height = 500,
            margin=dict(
                l=20,
                r=20,
                b=25,
                t=20,
                pad=5
            ),
        )
        outfig = [html.H3(value[0]),
                  dcc.Graph(
                        id = 'heatmap0',
                        figure = fig,
                        config = {
                            'toImageButtonOptions':{
                                'filename': 'degapp_heatmap_'+str(value[0])+'_alpha' + str(signif_cutoff),
                                'format': 'svg'
                            }
                        }
                 )]
    if debug == True:
        fig.show()
    return outfig

# ### Callback Gene
@app.callback(
    Output({'type':'gene-go-container', 'index': MATCH}, 'page_size'),
    [Input({'type':'select_page_size', 'index': MATCH}, 'value')])
def update_go_pagesize(value):
    return value

@app.callback(
    Output('gene-info-container', 'children'),
    Input('hm_gene_store', 'data'))
def update_info_panel(genes):
    if genes is None:
        return dbc.Row(no_selection_instruction)
    elif type(genes) == str:
        genes = [genes]
    
    if len(genes)==0:
        return dbc.Row(no_selection_instruction)
    else:
        qry_res = query_genes(genes)
        panel_ls = []
        for i in range(len(genes)):
            panel_ls.append(make_geneinfo_container(genes[i], i, qry_res= qry_res))
            panel_ls.append(html.P(style = {'backgroundColor':'#f9f9f9'}))

        return panel_ls


app.layout = [
    header.ui,
    header.help_page.ui,

    app_content
]

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8050, debug=True)
