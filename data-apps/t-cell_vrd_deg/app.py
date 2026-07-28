"""
app.py is a template for Allen Institute members to quickly, 
efficiently, and uniformly develop Dash apps
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "libs" / "allen_dash_modules" / "src"))

from dash import Dash

import allen_dash_modules as adm

from components.create_deg_components import create_deg_components
from ui.create_deg_layout import create_deg_layout
from callbacks.create_deg_callbacks import create_deg_callbacks

# parse arguments
parser = argparse.ArgumentParser(
    description=
"""
Allen Dash Modules (ADM) is a package that helps developers at the Allen Institute write Dash apps
more efficiently and generalizes their Dash apps to multiple datasets.
"""
)
parser.add_argument("-c", "--config", required=False, help="Path to an ADM-compatible JSON configuration file.")
parser.add_argument("-s", "--schema", required=False, help="Path to a JSON schema to validate configuration against.")
parser = parser.parse_args()

# run program
if __name__=='__main__':

    # create configuration singleton
    config = adm.ConfigSingleton(
        config_file="config/tcell-vrd.json",
        schema_file=parser.schema
    )

    # initiate dash app here
    app = Dash(
        name=config.dashboard.project_title,
        external_stylesheets=config.dashboard.stylesheets,
        update_title=""
    )
    app._favicon = "images/aifi_logo.png"

    # 1. create deg components
    deg_components = create_deg_components(config)

    # 2. define app layout with deg components
    app.layout = create_deg_layout(config, components=deg_components)

    # 3. define deg callbacks with deg components
    create_deg_callbacks(config, components=deg_components)

    app.run(host='0.0.0.0', port=8050, debug=False)
