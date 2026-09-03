"""
modal page to allow users to configure plot sizes
and miscellaneous components/behaviors by clicking the gear icon
"""

from os import PathLike

from dash_bootstrap_components import Modal, ModalBody

from ..utils import load_markdown

class StartupModal():
    """
    A class representing modal to give a brief overview of the project and orient users to the data app.

    Attributes
    ----------
    ui : dash_bootstrap_components.Modal
        The Modal component that gets placed into the app's layout.

    Methods
    -------
    build_startup_modal(markdown_file) : allow users to define the content of the startup modal

    """
    def __init__(self):
        self.ui:    Modal
        self.id:    str = 'startup-modal'

    def build_startup_modal(self, markdown_file: PathLike) -> Modal:
        """
        Constructs the startup modal using a custom markdown file. Works exactly as the HelpPanel class and
        the build_help_panel method. To get a title, write it in markdown using the headers.

        Parameters
        ----------
        markdown_file : PathLike
            Path to the Markdown file containing the modal text.

        Returns
        -------
        dbc.modal : list of dash components
            Same content as the input with extra pre-defined properties to adhere to ADM standards.
            App content can be saved into the layout for data visualization.
        """

        modal_body = ModalBody(load_markdown(markdown_file), className='markdown-body')

        modal = Modal(
            children = modal_body,
            id = self.id,
            is_open=True,
            size='lg',
            fade=True,
            centered=False
        )

        return modal
