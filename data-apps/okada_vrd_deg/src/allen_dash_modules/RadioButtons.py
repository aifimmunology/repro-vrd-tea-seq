from dash import dcc, html
import dash_bootstrap_components as dbc

class RadioButtons:

    def __init__(self, label: str, **kwargs):

        self.dd:    dcc.Dropdown
        self.ui:    dbc.InputGroupText
        self.id:    str = kwargs['id']

        # label
        radio_label = html.B(label, style={
                'fontSize': '1.1rem'
                })

        # radio button options
        radio_buttons = dbc.RadioItems(
            **kwargs,
            className="btn-group",
            inputClassName="btn-check",
            labelClassName="btn btn-outline-primary",
            labelCheckedClassName="active"
        )

        # everything together
        radio_group_children = dbc.InputGroupText(
                [
                    radio_label,
                    radio_buttons
                ]
            )

        self.radio_buttons = radio_buttons
        self.ui = radio_group_children
