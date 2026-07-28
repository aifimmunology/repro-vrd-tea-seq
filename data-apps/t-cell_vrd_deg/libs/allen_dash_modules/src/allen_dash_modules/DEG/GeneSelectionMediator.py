"""
This mediator is important to coordinate user click selection within the DEG volcano plot
along with the value of a dropdown. These two visual components can overwrite the value of
the dcc.Store "gene-selection-store", which can be used to affect additional downstream elements.
"""

from dash import callback, Input, Output, ctx, no_update

class GeneSelectionMediator:
    def __init__(self):
        pass

    def mediate(self, deg_plot, dropdown_gene):
        @callback(
            Output("gene-selection-store", 'data'),
            Output(dropdown_gene.id, "value"),
            Input(dropdown_gene.id, "value"),
            Input(deg_plot.id, 'clickData'),
            Input("ma_plot", "clickData"),
            Input("volcano_deg_grid", "cellClicked"),
            prevent_initial_call=True
        )
        def gene_selection_mediator(
            dropdown_gene_selection,
            click_figure,
            click_ma,
            click_table
        ):

            # update gene selection store via user click on table
            if 'grid' in ctx.triggered_id:
                if click_table is not None:
                    if click_table['colId'].lower() == 'gene':
                        return click_table['value'], click_table['value']
                    else:
                        return no_update, no_update
                else:
                    return no_update, no_update
                
            if ctx.triggered_id == 'ma_plot':
                selected_gene=click_ma['points'][0]['customdata'][0]
                return selected_gene, selected_gene

            # return None values if Gene dropdown is cleared
            if ctx.triggered_id == 'dropdown_gene':
                if dropdown_gene_selection is None:
                    return None, None
                else:
                    return dropdown_gene_selection, dropdown_gene_selection

            # if user clicks a point, then change the state of gene-selection-store
            if ctx.triggered_id == deg_plot.id:
                selected_gene=click_figure["points"][0]["customdata"][0]
                return selected_gene, selected_gene
