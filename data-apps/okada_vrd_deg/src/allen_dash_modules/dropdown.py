from dash import dcc, html
import dash_bootstrap_components as dbc

class DropDown:

    def __init__(self, id, label: str, **kwargs):

        self.dropdown:      dcc.Dropdown
        self.ui:            dbc.InputGroupText
        self.id:            str | dict = id

        # label on the left-hand side of the dropdown
        dropdown_label = html.B(
            label,
            style={
                'fontSize': '1rem',
                'marginRight': '0.5rem'
            })

        # dcc.Dropdown on the right-hand side
        dropdown = dcc.Dropdown(
            id=id,
            className='dbc',
            style={
                'width': '100%',
                'fontSize': '0.90rem',
                'textAlign': 'left'
            },
            **kwargs,
        )

        # assemble components via dbc.InputGroupText
        input_group_children = dbc.InputGroupText([
            dropdown_label,
            dropdown
        ])

        self.dropdown = dropdown
        self.ui = input_group_children
