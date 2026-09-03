"""
This module defines the HelpPanel class and a related function for displaying
an off-canvas sidebar with help content in a Dash app. The panel is toggled
on or off when the user clicks a help ("?") button.

The panel content is loaded from a Markdown file and rendered using Dash's
Markdown component, with optional MathJax support.
"""

from os import PathLike

from dash import Input, Output, State, callback
from dash_bootstrap_components import Offcanvas

from .utils import load_markdown

class HelpPanel:
    """
    A class representing an off-canvas help panel for displaying application help content.

    Attributes
    ----------
    ui : dash_bootstrap_components.Offcanvas
        The off-canvas UI component that holds the help content.

    Methods
    -------
    build_help_panel(markdown_file)

    """

    def __init__(self):
        """
        Initializes the HelpPanel object with default values.
        """
        self.ui: Offcanvas
        self.component_id: str = "ui-help"
        self.toggle_button_id: str = "toggle-help-page"

    def build_help_panel(self, markdown_file: PathLike):
        """
        Constructs the help panel UI element and assigns it to `self.ui`.

        Parameters
        ----------
        markdown_file : PathLike
            Path to the Markdown file containing the help text.

        Returns
        -------
        dbc.Offcanvas
            An Offcanvas page with text content written in Markdown and formatted
            by the CSS selectors 'offcanvas, offcanvas-help, markdown-body'

        """

        body = load_markdown(markdown_file)

        help_page = Offcanvas(
            children=body,
            id=self.component_id,
            is_open=False,
            scrollable=True,
            close_button=False,
            placement="end",
            className="offcanvas-help",
        )

        return help_page

    def help_callback(self, component_id, component_button):
        """
        Registers a Dash callback to toggle the visibility of the help panel.

        Parameters
        ----------
        component_id : str
            The ID of the Offcanvas page that is found within this class. 'ui-help'

        component_button : str
            The ID of the button that triggers found in this class. 'toggle-help-page'

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
