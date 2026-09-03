"""
This module defines the AppHeader class, which constructs a navigation header
for a Dash app using Dash Bootstrap Components. The header includes the app title,
external links (e.g., institute and study URLs), optional settings and help buttons,
and a light/dark mode toggle switch.
"""

from dash_iconify import DashIconify
from dash import Input, Output, html, clientside_callback
from dash_bootstrap_components import Container, Navbar, NavbarBrand, Switch, NavLink, Row, Col

from .app_settings import SettingsPanel
from .app_help import HelpPanel

def header_button(id, icon, height=22):
    """
    Helper function to create a styled icon button using DashIconify.

    Parameters
    ----------
    id : str
        The ID to assign to the button.
    icon : str
        The icon name for DashIconify.

    Returns
    -------
    dash.html.Button
        A styled button with the specified icon.
    """
    button = html.Button(
        id=id,
        children=[DashIconify(icon=icon, height=height, style={"color": "white"})],
        style={"background": "transparent", "border": "none", "cursor": "pointer"},
    )

    return button


class AppHeader:
    def __init__(
        self,
        toggle_theme: bool,
        toggle_settings: bool,
        toggle_help: bool
    ):
        """
        Attributes
        ----------

        toggle_help : bool
            Include a help button, offCanvas page, and associating callback to open the page.

        toggle_settings : bool
            Include a settings button, offCanvas page, and associating callback to open the page.

        toggle_theme : bool
            Include a light and dark control into the app header with component id of "toggle-dark-mode".

        Notes
        -----
        A class for creating a customizable header bar for Dash applications. Additional
        arguments can be used to configure the look and behavior of the navigation bar.
        Component styling can be adjusted using CSS. \n

        This object uses class composition techniques to build the header and is composed of a
        help page and settings page. To access and build an offCanvas component for these pages,
        access them using the following method:\n

        create a help page\n
        header.help_page.build_help_panel(
            component_id="help-page",
            markdown_file=config.dashboard.help_page
        )

        create a settings page\n
        my_addition_components = html.Div([dbc.Slider(...), dbc.Button])\n
        header.settings_page.extend(my_addition_components)\n

        """
        self.ui: Navbar
        self.component_id: str = "ui-header"

        self.toggle_theme: bool = toggle_theme
        self.theme_switch_id: str = "toggle-dark-mode"

        self.toggle_settings: bool = toggle_settings
        self.settings_button_id: str = "toggle-settings-page"
        self.settings_page: SettingsPanel

        self.toggle_help: bool = toggle_help
        self.help_button_id: str = "toggle-help-page"
        self.help_page: HelpPanel

        if self.toggle_theme:
            self.theme_clientside_callback(self.theme_switch_id)

        if self.toggle_help:
            self.help_page = HelpPanel()
            self.help_page.help_callback(
                component_id=self.help_page.component_id,
                component_button=self.help_page.toggle_button_id,
            )

        if self.toggle_settings:
            self.settings_page = SettingsPanel()
            self.settings_page.settings_callback(
                component_id=self.settings_page.component_id,
                component_button=self.settings_page.toggle_button_id,
            )

    def build_app_header(self, title: str, institute_url: str, study_url: str):
        """
        Constructs and stores the app header layout in `self.ui`.

        Parameters
        ----------
        title : str
            The title displayed in the app header.
        institute_url : str
            URL linking to the developer's organization or institution.
        study_url : str
            URL linking to the study or project documentation.
        id : str, optional
            ID for the top-level Navbar element (default is 'ui-header').
        self.settings_page : bool, optional
            Whether to include a settings icon button (default is False).
        settings_button_id : str, optional
            ID for the settings button (default is 'button-header-settings').
        toggle_theme : bool, optional
            Whether to include a dark/light mode toggle (default is True).
        toggle_theme_id : str, optional
            ID for the toggle switch (default is 'toggle-dark-mode').
        help_button : bool, optional
            Whether to include a help icon button (default is True).
        help_button_id : str, optional
            ID for the help button (default is 'button-header-help').

        Returns
        -------
        None
        """

        if self.toggle_settings:
            settings_obj = [
                " | ",
                header_button(self.settings_button_id, "carbon:settings"),
            ]
        else:
            settings_obj = []

        if self.toggle_theme:
            mode_obj = [
                " | ",
                html.Div(
                    [
                        header_button("moon", "ph:moon"),
                        Switch(
                            id=self.theme_switch_id,
                            value=True,
                            input_style={
                                "width": "3rem",
                                "height": "1.5rem",
                            },
                            persistence=True,
                            persistence_type='local'
                        ),
                        header_button("sun", "ph:sun"),
                    ],
                    style={"display": "inline-flex"},
                ),
            ]
        else:
            mode_obj = []

        if self.toggle_help:
            help_obj = [
                " | ",
                header_button(self.help_button_id, "material-symbols:help-outline"),
            ]
        else:
            help_obj = []

        header = Navbar(
            id=self.component_id,
            children=[
                Container(
                    fluid=True,
                    children=[
                        Row(
                            className="w-100 align-items-center g-0",
                            justify="between",
                            children=[

                                # LEFT GROUP: logo, "allen institute /", and app title (left-aligned, middle-aligned)
                                Col(
                                    width=9,
                                    className="d-flex align-items-center justify-content-start",
                                    children=[

                                        # link with logo
                                        NavLink(
                                            href=institute_url,
                                            target="_blank",
                                            external_link=True,
                                            className="p-1 d-flex align-items-center",
                                            children=[
                                                html.Img(
                                                    src="assets/images/AI_lens_white.svg",
                                                    height=36,
                                                    className="me-1",
                                                ),
                                            ],
                                        ),

                                        NavLink(
                                            href=institute_url,
                                            target="_blank",
                                            external_link=True,
                                            className="p-1 d-flex align-items-center",
                                            children=html.Span("allen institute /", style={"fontSize": 20}),
                                        ),

                                        NavLink(
                                            href=institute_url,
                                            target="_blank",
                                            external_link=True,
                                            className="p-1 d-flex align-items-center",
                                            children=html.Span("immunology /", style={"fontSize": 20}),
                                        ),

                                        # project title
                                        NavLink(
                                            href=study_url,
                                            target="_blank",
                                            external_link=True,
                                            className="p-1 d-flex align-items-center",
                                            children=html.Span(title, style={'fontSize': 20})
                                        )
                                    ],
                                ),

                                # RIGHT GROUP: study overview, settings, theme toggle, help (right-aligned)
                                Col(
                                    width=3,
                                    className="d-flex align-items-center justify-content-end gap-1",
                                    children=[
                                        html.A(
                                            href=study_url,
                                            target="_blank",
                                            children=[
                                                header_button(
                                                    "beaker",
                                                    "fluent:document-one-page-beaker-16-regular",
                                                ),
                                                " Study Overview",
                                            ],
                                            style={
                                                "textDecoration": "none",
                                                "color": "white",
                                                "marginRight": "5px",
                                            },
                                        )
                                    ]
                                    + settings_obj
                                    + mode_obj
                                    + help_obj,
                                ),
                            ],
                        ),
                    ],
                    className="navbar-container",
                )
            ],
            className="navbar",
            color="white",
        )

        return header

    def theme_clientside_callback(self, toggle):
        """
        Sets up a clientside callback to toggle between light and dark mode
        using the specified Switch component.

        Parameters
        ----------
        toggle : str
            The ID of the Switch component that controls theme switching.

        Returns
        -------
        None
        """
        clientside_callback(
            """
            (switchOn) => {
            document.documentElement.setAttribute("data-bs-theme", switchOn ? "light" : "dark");
            return window.dash_clientside.no_update
            }
            """,
            Output(toggle, "id"),
            Input(toggle, "value"),
        )
