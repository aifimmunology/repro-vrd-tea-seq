from allen_dash_modules import AppHeader
import dash_bootstrap_components
import dash

title = "Test Title"
aifi_url = "https://alleninstitute.org/division/immunology/"
study_url = "https://explore.allenimmunology.org/"

def test_init_header():
    header = AppHeader()

    assert header.ui is None

def basic_header_tests(header):

    assert isinstance(header.ui, dash_bootstrap_components.Navbar)
    # Navbar.Container Should have 3 objects: Logos, Title, and Study Link/Settings
    assert len(header.ui.children[0].children) == 3
    # Navbar.Container[0].Div[0].A[0] Should point to the Allen Institute URL
    assert header.ui.children[0].children[0].children[0].href == aifi_url
    # Navbar.Container.Div[1] Should have our title
    assert header.ui.children[0].children[1].children[0].children == title
    # Navbar.Container.Div[2].A[0] should point to the Study URL
    assert header.ui.children[0].children[2].children[0].href == study_url

def test_basic_header():
    header = AppHeader()

    header.build_app_header(
        title = title,
        study_url = study_url
    )

    basic_header_tests(header)  
   
def test_help_header():
    header = AppHeader()

    header.build_app_header(
        title = title,
        study_url = study_url,
        help_button = True
    )
     
    # Ensure basic tests still pass
    basic_header_tests(header)

    # Navbar.Container.Div[2].Div[2] should be a button
    assert isinstance(header.ui.children[0].children[2].children[2], dash.html.Button)
    # Navbar.Container.Div[2].Div[2] The button's id should be the default value
    assert header.ui.children[0].children[2].children[2].id == 'button-header-help'

def test_help_custom_header():
    header = AppHeader()

    custom_id = 'custom_help_id'

    header.build_app_header(
        title = title,
        study_url = study_url,
        help_button = True,
        help_button_id = custom_id
    )
     
    # Ensure basic tests still pass
    basic_header_tests(header)

    # Navbar.Container.Div[2].Div[2] should be a button
    assert isinstance(header.ui.children[0].children[2].children[2], dash.html.Button)
    # Navbar.Container.Div[2].Div[2] The button's id should be the custom value
    assert header.ui.children[0].children[2].children[2].id == custom_id

def test_settings_header():
    header = AppHeader()

    header.build_app_header(
        title = title,
        study_url = study_url,
        settings_button = True
    )
     
    # Ensure basic tests still pass
    basic_header_tests(header)

    # Navbar.Container.Div[2].Div[2] should be a button
    assert isinstance(header.ui.children[0].children[2].children[2], dash.html.Button)
    # Navbar.Container.Div[2].Div[2] The button's id should be the default value
    assert header.ui.children[0].children[2].children[2].id == 'button-header-settings'

def test_settings_custom_header():
    header = AppHeader()

    custom_id = 'custom_settings_id'

    header.build_app_header(
        title = title,
        study_url = study_url,
        settings_button = True,
        settings_button_id = custom_id
    )
     
    # Ensure basic tests still pass
    basic_header_tests(header)

    # Navbar.Container.Div[2].Div[2] should be a button
    assert isinstance(header.ui.children[0].children[2].children[2], dash.html.Button)
    # Navbar.Container.Div[2].Div[2] The button's id should be the custom value
    assert header.ui.children[0].children[2].children[2].id == custom_id
