# Allen-Dash-Modules

1. [What is this project about?](#what-is-this-project-about)

    a. [Application 1: Develop Dash apps quickly and efficiently](#application-1-develop-dash-apps-quickly-and-efficiently)

    b. [Application 2: Apply a Dash app across multiple datasets](#application-2-apply-a-dash-app-across-multiple-datasets)

2. [Frequently Asked Questions](#frequently-asked-questions)

3. [Quick Start](#quick-start)

## What is this project about?

Allen Dash Modules (ADMs) is a collection of modules to **standardize visual and behavioral components of Dash apps for Allen Institute members to quickly and efficiently write data apps**. By encapsulating frequently used components, plots, and functions into Python classes, ideally users can pull and adapt re-usable components to fit a wide variety of applications in bioinformatics and data science. Another goal of this package is to **enable developers to generalize one data app across multiple datasets**. This approach can reduce redundant efforts to code, enhance reproducibility of results, and allow developers to focus on component interactions and behaviors. The implementations of the two goals above are further expanded below:

### Application 1: Develop Dash apps quickly and efficiently

We packaged a variety of plots and callbacks into classes. For example, the following code snippet allows developers to create a standardized header that fits the Allen Institute themes, while allowing configuration of specific parameters.

```python
from dash import Dash, html

import dash_bootstrap_components as dbc

import allen_dash_modules as adm

# create header object
header = adm.AppHeader(
    toggle_help=True,
    toggle_settings=True,
    toggle_theme=True
)

# build app header
header.ui = header.build_app_header(
    title="Test data app",
    institute_url="https://alleninstitute.org/division/immunology/",
    study_url="https://explore.allenimmunology.org/"
)

# build help page by providing a markdown file
header.help_page.ui = header.help_page.build_help_panel(
    markdown_file="assets/markdown/template_help_page.md"
)

# build settings page by building an Offcanvas page
header.settings_page.ui = header.settings_page.build_settings_panel(
    children=dbc.Offcanvas(
        [
            html.H2("Dashboard Configuration"),
            html.Span("Include additional components below!")
        ]
    )
)

# create a startup modal using a markdown file
startup_modal = adm.StartupModal()
startup_modal.ui = startup_modal.build_startup_modal(
    markdown_file="assets/markdown/template_startup_modal.md"
)

# establish stylesheets
stylesheets = [
    "https://cdn.jsdelivr.net/npm/bootstrap@5.3.6/dist/css/bootstrap.min.css",
    "https://storage.googleapis.com/aifi-static-assets/hise-style.css",
    "assets/css/dbc.min.css",
    "assets/css/sidebar.css",
    "assets/css/markdown.css",
    "assets/css/startup_modal.css",
    "assets/css/layout.css"
]

# build dash app
app = Dash(external_stylesheets=stylesheets, update_title="")
app._favicon = "images/aifi_logo.png"

# create a mock layout with a sidebar
mock_sidebar = html.Div(
    [
        "sidebar"
    ],
    style={
        'width': '25rem',
        'height': '90vh',
        'padding': '1rem',
        'backgroundColor': 'rgba(143, 143, 143, 0.5)',
        'marginRight': '1rem'
    }
)

# mock layout with main content area
mock_main_content = html.Div(
    [
        "main content"
    ],
    style={
        'width': '100%',
        'height': '90vh',
        'padding': '1rem',
        'backgroundColor': 'rgba(0, 102, 147, 0.5)',
    }
)

# assemble sidebar and main content into one Div with custom 'container' class
app_content = html.Div(
    [
        mock_sidebar,
        mock_main_content
    ],
    className='container'
)

# assemble pre-made components into the layout
app.layout = [
    header.ui,
    header.settings_page.ui,
    header.help_page.ui,
    startup_modal.ui,

    app_content
]

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8050, debug=True)

```

The above code would create a navigation bar on the top with two offCanvas pages to store settings configuration and a help page shown below:

![alt text](assets/images/template.png)

To add additional features to your app, several commonly used visual components have been pre-packaged to quickly develop data apps in a standardized way as displayed below:

| Class                    | Folder     | Description                                                                                                                                                                            |
|--------------------------|------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| StartupModal             | components | Orient users on how to use your app and what goals can be accomplishesd with this tool via a Markdown file.                                                                            |
| ConfigSingleton          | components | Save a JSON configuration file as a Python class and dictionary in one object. Allows devs to access specific values using dot notation.                                               |
| GenePanel                | components | GenePanel class to build a Dash component for displaying gene information retrieved from the MyGene.info API.                                                                          |
| RadioButtons             | components | Make multiple buttons next to each other.                                                                                                                                              |
| DropDown                 | components | Make a stylized dropdown menu that responds to light and dark mode.                                                                                                                    |
| AppHeader                | header     | Make a standardized navigation bar at the top of the screen with optional pages for settings, help message written via Markdown file, and light/dark switch toggle.                    |
| VolcanoPlot              | DEG        | Display up and down-regulated differential genes with respect to log fold change and p-value.                                                                                          |
| MA Plot                  | DEG        | Display up and down-regulated differential genes with respect to log fold change and means of gene expression between two conditions.                                                  |
| DEG                      | DEG        | Holds all the components in a DEG app for quick retrieval to build the app layout using techniques of [class composition](https://en.wikipedia.org/wiki/Composition_over_inheritance). |
| GeneSelectionMediator    | DEG        | Listen to several input sources and eventually write to a dcc.Store that holds either a selected gene (str) or a None value.                                                           |
| PathwaySelectionMediator | DEG        | Listen to several input sources and eventually write to a dcc.Store that holds either a selected pathway (str) or a None value.                                                        |
| BoxPlot                  | plots      | Create a boxplot to optionally show perturbation status in paired samples. Used for differential gene expression.                                                                      |
| MAPlot                   | plots      | Create an MA plot to show differential gene expression via log fold change and mean gene expression values.                                                                            |
| Volcano                  | plots      | Create a Volcano plot to show differential gene expression via log fold change and p-value.                                                                                            |

### Application 2: Apply a Dash app across multiple datasets

After developing a satisfactory app, consider applying your app to a different dataset to further increase its utility! This can be facilitated by using the pre-installed configuration file `config/template.json`:

```json
{
    "dashboard": {
        "institute_homepage": "https://alleninstitute.org/division/immunology",
        "project_homepage": "https://explore.allenimmunology.org/",
        "project_title": "Test data app",
        "help_page_content": "assets/markdown/template_help_page.md",
        "favicon": "assets/images/aifi_logo.png",
        "stylesheets": [
            "https://storage.googleapis.com/aifi-static-assets/hise-style.css",
            "https://cdn.jsdelivr.net/npm/bootstrap@5.3.6/dist/css/bootstrap.min.css",
            "assets/css/dbc.min.css",
            "assets/css/sidebar.css",
            "assets/css/markdown.css",
            "assets/css/startup_modal.css",
            "assets/css/layout.css"
        ],

        "schema": "tests/schema.json",

        "startup_modal": true,
        "startup_modal_content": "assets/markdown/template_startup_modal.md"

    }
}
```

The main advantages of using an input JSON file is to configure paths to a data source, adjust column names to accommodate projects, and optionally include components or modify behaviors. Input JSON files can be automatically checked for valid data types using tools like [jsonschema](https://python-jsonschema.readthedocs.io/en/stable/). This approach should allow developers to focus on component interactions and programming logic, which can be applied to different projects. ADM have been previously applied to a scRNA-Seq dataset of human blood samples perturbed by 90 cytokines generated by Parse Biosciences and the configuration file is shown below. The entire codebase for that project can be found [here](https://github.com/aifimmunology/allen-dash-modules/tree/2b3ee87eb851cb66e467bb2c69f7eced05b2a17c).

<details>
    <summary>Click to expand</summary>

```markdown

{
    "dashboard": {
        "institute_homepage": "https://alleninstitute.org/division/immunology",
        "project_homepage": "https://apps.allenimmunology.org/aifi/resources/parse-10m-cytokines",
        "project_title": "Perturbation of 10M PBMCs with 90 Cytokines with Parse Biosciences GigaLab",
        "help_page": "assets/help_page.md",

        "schema": "tests/schema.json",

        "gene_panel": true,

        "startup_modal": true,
        "startup_modal_title": "Allen Institute for Immunology | Differential Gene Expression Explorer",
        "startup_modal_content": "assets/startup_modal.md"

    },

    "experiment": {
        "covariates": ["Cytokine", "Cell Type", "Censoring"],
        "perturbation": true,
        "control": "PBS"
    },

    "data": {

        "deg": {
            "specifications": {
                "data_type": "deg",
                "file_type": "parquet",

                "ma_plot": true,
                "deg_by_covariates": true,
                "summaries_plot": true,

                "p_col": "padj",
                "es_col": "log2fc",
                "means_col": "means",
                "store": "volcano-deg-df"
           },
            "files": {
                "adjp": "data/deg/adjp.parquet",
                "log2fc": "data/deg/log2fc.parquet",
                "meta": "data/deg/meta.parquet",
                "means": "data/deg/mean.parquet",
                "summaries": "data/deg/summaries.parquet"
            }
        },

        "Hallmark": {
            "specifications": {
                "data_type": "gsea",
                "file_type": "parquet",
                "p_col": "adjP",
                "es_col": "NES",
                "store": "volcano-gsea-df"
            },
            "files": {
                "results": "data/gsea/hallmark.parquet",
                "meta": "data/gsea/hallmark_meta.parquet",
                "summaries": "data/gsea/hallmark_summaries.parquet"
            }
        },

        "Reactome": {
            "specifications": {
                "data_type": "gsea",
                "file_type": "parquet",
                "p_col": "adjP",
                "es_col": "NES",
                "store": "volcano-gsea-df"
            },
            "files": {
                "results": "data/gsea/reactome.parquet",
                "meta": "data/gsea/reactome_meta.parquet",
                "summaries": "data/gsea/reactome_summaries.parquet"
            }
        },

        "pseudobulk": {
            "specifications": {
                "data_type": "pseudobulk",
                "file_type": "parquet",
                "bulk_label": "Donor",
                "store": "pb-df"
            },
            "files": {
                "genes": "data/pseudobulk/genes.parquet",
                "meta": "data/pseudobulk/meta.parquet",
                "norm": "data/pseudobulk/norm.parquet"
            }
        }
    }
}

```

</details>

## Quick Start
In order to run the `template.py` file, ensure that either [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install) or the Python [uv](https://docs.astral.sh/uv/guides/install-python/) package managers are installed. Then input only one set of the following commands into a terminal:

```bash
# create anaconda environment
$ mamba env create -n adm-template -f environment.yml

# or install packages using uv
$ uv sync
```

Then run your application using the following command in the command line:

```bash
# clone the repo
$ git clone git@github.com:aifimmunology/allen-dash-modules.git .

# move into the directory
$ cd allen-dash-modules

# move the script to the top-level directory!
$ mv apps/template.py .

# if running miniconda, run the script like this
$ python template.py

# if using uv, run the script like this
$ uv run template.py
```

To personalize your program, feel free to rename template.py and add additional components and callbacks into the file. Another option would be to create a new function in a separate file and then call that function. Please stay tuned for future updates on the best practices on how to use this program.

## Frequently asked questions

### Why are my popovers/tooltips not behaving consistently?

The popover/tooltip must object must be in the same list as your target object.

### Why is my app not stylized at all?

Make sure the script being run is located at the top-level of the working directory (above the assets/ folder). This can mean moving the template.py file from the apps/ folder up one directory.

## Resources

Meijer, P., Howard, N., Liang, J., Kelsey, A., Subramanian, S., Johnson, E., Mariz, P., Harvey, J., Ambrose, M., Tereshchenko, V., et al. (2025). Provide proactive reproducible analysis transparency with every publication. R. Soc. Open Sci. 12, 241936. https://royalsocietypublishing.org/doi/10.1098/rsos.241936