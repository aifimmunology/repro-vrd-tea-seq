from dash import Input, Output, State, callback, Dash, ALL

from .GeneSelectionMediator import GeneSelectionMediator
from .PathwaySelectionMediator import PathwaySelectionMediator
from ..components import ConfigSingleton

class DEG:

    def __init__(self, config: ConfigSingleton, components: dict):

        # save each component using their component id.
        for component_id, component in components.items():
            if isinstance(component, tuple):
                for sub_component in component:
                    # use the sub-component's id attribute as the name
                    if isinstance(sub_component.id, str):
                        setattr(self, sub_component.id, sub_component)
                    # do the same thing but for id's that fit the pattern-matching callback format 
                    elif isinstance(sub_component.id, dict):
                        setattr(self, sub_component.id["index"], sub_component)
            else:
                setattr(self, component_id, component)

        # establish mediator classes to define app behavior

        # gene selection mediator. this is a required mediator!
        gene_mediator = GeneSelectionMediator()
        gene_mediator.mediate(
            deg_plot=self.volcano_deg,
            dropdown_gene=self.dropdown_gene
        )

        # optional gsea mediator
        gsea_sources = [ i for i in config.dict["data"].keys() if config.dict["data"][i]["specifications"]["data_type"] == "gsea" ]
        if gsea_sources:
            pathway_mediator = PathwaySelectionMediator()
            pathway_mediator.mediate(
                gsea_plot=self.volcano_gsea,
                dropdown_covariates=components["dropdown_covariates"],
                dropdown_pathway=self.dropdown_pathway
            )

        # optional mediator for volcano plot of 1 gene across all covariates
        if "volcano_by_covariates_plot" in dir(self):
            @callback(
                Output({"type": "covariate", "index": ALL}, "value", allow_duplicate=True),
                State("covariate-values", "data"),
                Input("volcano_by_covariates_plot", "clickData"),
                prevent_initial_call=True
            )
            def VolcanoWrtCovariatesMediator(values, click):
                if click:
                    covariates = click["points"][0]["customdata"][0].split(": ")
                    covariates.extend([values[-1]])
                    return covariates

        if "deg_summaries" in dir(self):
            @callback(
                Output({"type": "covariate", "index": ALL}, "value", allow_duplicate=True),
                State("covariate-values", "data"),
                Input("deg_summaries", "clickData"),
                prevent_initial_call=True
            )
            def DegSummariesCovariatesMediator(values, click):
                if click:
                    values[1] = click['points'][0]['label']
                    return values
        # TODO: mitigate manual indexing by setting the covariates as a dictionary where key is label and value is value
        # NOTE: current indexing scheme fits the Parse 10M dataset
