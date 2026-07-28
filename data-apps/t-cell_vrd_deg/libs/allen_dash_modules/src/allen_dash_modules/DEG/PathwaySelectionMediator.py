"""
This mediator is important to coordinate user click selection within the GSEA volcano plot
along with the value of a dropdown. These two visual components can overwrite the value of
the dcc.Store "pathway-selection-store", which further changes the highlighted points in the
DEG volcano plot.
"""

from base64 import b64decode
from json import loads
from dash import callback, Input, Output, State, ctx, no_update, ALL

from ..components import DropDown

def reverse_lookup(dictionary: dict[str, str], target_value):
    for key, value in dictionary.items():
        if value == target_value:
            return key

class PathwaySelectionMediator:
    def __init__(self):
        pass

    def mediate(self,
        gsea_plot_id,
        gsea_grid_id,
        gsva_plot_id,
        gsva_grid_id,
        dropdown_covariates:    list[DropDown],
        dropdown_pathway:       DropDown
    ):
        
        plot_ids = [gsea_plot_id, gsva_plot_id]
        table_ids = [gsea_grid_id, gsva_grid_id]

        @callback(
            Output("pathway-selection-store", 'data'),
            Output(dropdown_pathway.id, 'value'),

            Input(gsea_plot_id, 'clickData'),
            Input(gsea_grid_id, "cellClicked"),

            Input(gsva_plot_id, 'clickData'),
            Input(gsva_grid_id, "cellClicked"),

            Input(dropdown_pathway.id, 'value'),
            State("pathway-label-store", 'data'),
            State("sticky-pathway-checkbox", "value"),

            [Input(dd.id, 'value') for dd in dropdown_covariates],
            prevent_initial_call=True
        )
        def pathway_selection_mediator(
            gsea_click_data,
            gsea_click_table,
            gsva_click_data,
            gsva_click_table,
            dropdown_pathway,
            pathway_label_store,
            sticky_pathway_checkbox,
            click_table,
            *args
        ):

            # this command flow could be a series of functions if you want,
            # but should be understandable as-is.

            # decode serialized dictionary
            pathway_label_store = loads(b64decode(pathway_label_store.encode('utf-8')).decode('utf-8'))

            # define covariate-related dd's AKA all but the last one
            covariate_dd_ids = tuple(dd.id for dd in dropdown_covariates)

            # if change in covariate dd, reset pathway selection
            # optionally, make the pathway sticky across selections with a switch.
            if sticky_pathway_checkbox:
                if ctx.triggered_id in covariate_dd_ids:
                    if ctx.inputs["dropdown_pathway.value"] is not None:
                        pathway_label=ctx.inputs["dropdown_pathway.value"]
                        pathway_code=reverse_lookup(pathway_label_store, ctx.inputs["dropdown_pathway.value"])
                        return pathway_code, pathway_label
                    else:
                        return None, None
            else:
                if ctx.triggered_id in covariate_dd_ids:
                    return None, None

            if ctx.triggered_id in table_ids:
                if ctx.triggered_id == gsea_grid_id:
                    click_table = gsea_click_table
                elif ctx.triggered_id == gsva_grid_id:
                    click_table = gsva_click_table
                else:
                    return no_update, no_update

                if click_table is not None:
                    if click_table['colId'].lower() == 'pathway label':
                        pathway_code=reverse_lookup(pathway_label_store, click_table['value'])
                        return pathway_code, click_table['value']
                    else:
                        return no_update, no_update
                else:
                    return no_update, no_update

            # Note: pathway_code is a computer-friendly identifier for pulling data in the back-end
            # Note: while pathway_label is a user-friendly string that is displayed on the UI.
            # Note: the following chunks uses the pathway_selection_store to direct results to each output.
                # {"PATHWAY_CODE": "Pathway Label"}
                # {"HALLMARK_APOPTOSIS": "APOPTOSIS"}

            # reset pathway selection store if pathway dd is cleared
            # otherwise if user makes a selection in pathway dd, set state
            if ctx.triggered_id == 'dropdown_pathway':
                if ctx.inputs["dropdown_pathway.value"] is None:
                    return None, None
                else:
                    pathway_label=ctx.inputs["dropdown_pathway.value"]
                    pathway_code=reverse_lookup(pathway_label_store, pathway_label)
                    return pathway_code, pathway_label

            # if user clicks a point, then change state and update pathway dropdown dd
            if ctx.triggered_id in plot_ids:
                if ctx.triggered_id == gsea_plot_id:
                    click_data = gsea_click_data
                elif ctx.triggered_id == gsva_plot_id:
                    click_data = gsva_click_data
                else:
                    return no_update, no_update
                
                pathway_label=click_data["points"][0]["customdata"][0]
                pathway_code=reverse_lookup(pathway_label_store, pathway_label)
                return pathway_code, pathway_label
