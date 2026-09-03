"""
Offcanvas page to allow users to configure plot sizes
and miscellaneous components/behaviors by clicking the gear icon
"""

from dash import Input, Output, State, callback
from dash_bootstrap_components import Offcanvas


class SettingsPanel:
    """
    A class representing an off-canvas panel to hold dbc/dcc components for configuration purposes.

    Attributes
    ----------
    ui : dash_bootstrap_components.Offcanvas
        The off-canvas UI component that holds the help content.

    Methods
    -------
    build_settings_panel(children) : allow users to define dash components inside a dbc.Offcanvas page

    """

    def __init__(self):
        self.ui: Offcanvas
        self.component_id: str = "ui-setting"
        self.toggle_button_id: str = "toggle-settings-page"

    def build_settings_panel(self, children: Offcanvas):
        """
        Constructs the settings panel Offcanvas component.

        Parameters
        ----------
        children : list of dash components
            The contents to be rendered inside the settings panel (e.g., sliders, dropdowns, divs, spans, etc.).

        Returns
        -------
        dbc.Offcanvas : list of dash components
            Same content as the input with extra pre-defined properties to adhere to ADM standards.
            App content can be saved into the layout for data visualization.
        """

        # make input Offcanvas properties explicit
        children.id = self.component_id
        children.scrollable = True
        children.close_button = False
        children.placement = "end"
        children.className = "offcanvas-settings"

        return children

    def settings_callback(self, component_id, component_button):
        """
        Registers a Dash callback to toggle the visibility of the settings panel.

        Parameters
        ----------
        toggle_button_id : str
            The ID of the button that toggles the settings panel.

        Returns
        -------
        None
        """

        @callback(
            Output(component_id, "is_open"),
            Input(component_button, "n_clicks"),
            [State(component_id, "is_open")],
        )
        def toggle_panel(n_clicks: int, is_open: bool):
            if n_clicks:
                return not is_open
            return is_open
