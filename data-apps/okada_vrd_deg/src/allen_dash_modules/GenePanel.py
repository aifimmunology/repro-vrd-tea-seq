"""
GenePanel class to build a Dash component for displaying gene information 
retrieved from the MyGene.info API.
Source of inspiration: https://github.com/aifimmunology/data-apps-vis/blob/development/apps/parse-10m-cyto-deg/components/gene_panel.py
"""
from mygene import MyGeneInfo
from dash_bootstrap_components import Card, CardHeader, CardBody, ListGroup, ListGroupItem
from dash import Input, Output, callback, html

class GenePanel():
    """
    GenePanel is a UI component that displays Gene 
    information pulled from the MyGene.info database.

    Attributes:
        id (str): Component ID for the gene selection panel.
        gene_info_id (str): Component ID for the gene info display container.
        feature (str): Label or type of feature (default: "Gene").
        ui (Component): Cached Dash layout component representing the panel.
    """
    def __init__(
        self,
        id:             str,
        gene_info_id:   str,
        output_div:     str = "gene_info_panel",
        input_store:    str = "gene-selection-store",
        feature:        str = "Gene"
    ):
        """
        Initialize the GenePanel object with layout and identifiers.

        Args:
            id (str): Unique ID for the gene panel container.
            gene_info_id (str): ID for the dynamic gene info content.
            feature (str): Type of entity being displayed (default: "Gene").
        """
        self.id = id
        self.gene_info_id = gene_info_id
        self.output_div = output_div
        self.input_store = input_store
        self.feature = feature

        self.ui = None

        self.build_gene_panel()

        @callback(
            Output(self.output_div, "children"),
            Input(self.input_store, "data")
        )
        def gene_panel_display(gene):

            if gene is None:
                info_panel=self.build_info_panel(gene="PDE4B")
            else:
                info_panel=self.build_info_panel(gene)
            return info_panel

    def build_gene_panel(self):
        """
        Builds the top-level layout (a Bootstrap card) for gene selection.

        This is a static wrapper UI that provides a container where gene 
        information can be dynamically rendered in `self.gene_info_id`.

        Sets:
            self.ui (Card): Dash component containing the gene info layout.
        """
        ui = Card([
            CardBody(id=self.gene_info_id)
        ], style={
            'height': '100%',
            'width': '100%'
        })

        self.ui = ui

    def query_mygene(
        self,
        gene: str,
        fields: list[str] = ['name','alias','summary']
    ):
        '''Perform mygene query with specific fields
        
        Params
        ------
        value list
            List of genes (gene symbols) to query from mygene
        myfields list
            List of field names to query. default ['name','alias','summary']
        
        '''
        mg = MyGeneInfo()

        mygene_res =  mg.querymany(
            qterms = gene,
            scopes = 'symbol',
            fields = fields,
            species = 'human',
            verbose=False
        )
        # why does the above function generate the following Warning?
        # Input sequence provided is already in string format. No operation performed
        
        mygene_res = mygene_res[0]

        # return None if no results from mygene
        if 'notfound' in mygene_res.keys():
            mygene_res = None
        else:

            # format incomplete fields, namely...

            # 1) no gene alias
            try:
                mygene_res["alias"]
            except KeyError:
                mygene_res['alias']="No gene aliases according to the MyGene database."

            # 2) only 1 alias available
            if isinstance(mygene_res["alias"], str) and len(mygene_res["alias"]) == 1:
                pass
            elif isinstance(mygene_res["alias"], list):
                mygene_res['alias'] = '; '.join(mygene_res['alias'])

            # 3) no summary results
            try:
                mygene_res["summary"]
            except KeyError:
                mygene_res["summary"] = "Summary unavailable according to the the MyGene database."

            return mygene_res

    def build_info_panel(self, gene: str):
        """
        Constructs a Bootstrap-styled info panel displaying gene metadata.

        Args:
            gene (str): Gene symbol used for querying MyGene.info.

        Returns:
            ListGroup: A Dash Bootstrap list of items showing gene metadata
                       or a warning if the gene is not found.
        """
        mygene_res = self.query_mygene(gene = gene)

        if mygene_res is None:
            ui = Card([
                    CardHeader(f"Gene Info for {gene}"),
                        ListGroup([
                            ListGroupItem(html.B(f'Gene Info for {gene}')),
                            ListGroupItem('Gene symbol not found in MyGene Database')
                        ])
                    ])
        else:
            ui = Card([
                    CardHeader(f"Gene Info for {gene}"),
                        ListGroup([
                            ListGroupItem([
                                html.B('Gene Name: '),
                                html.Span(mygene_res['name'])
                            ]),
                            ListGroupItem([
                                html.B('Gene Aliases: '),
                                html.Span(mygene_res['alias'])
                            ]),
                            ListGroupItem([
                                html.B('Gene Summary: '),
                                html.Span(mygene_res['summary'])
                            ])
                    ])
            ])

        return ui
