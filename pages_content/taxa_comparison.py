import streamlit as st
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.colors as pc
import plotly.graph_objects as go
import ast
import random
import numpy as np
import math

## Load data
@st.cache_data
def load_data():
    # Load data
    region_summary = pd.read_csv("data/region_summary.csv")
    virgo2_inventory = pd.read_csv('data/MAG_inventory_VIRGO2_021623_30Jul2024.txt', sep='\t')
    cluster_blast = pd.read_csv('data/cluster_blast.csv.gz')
    taxa_colors = pd.read_csv('data/VIRGO2_taxaKey.csv')
    return region_summary, virgo2_inventory, cluster_blast, taxa_colors

region_summary, virgo2_inventory, cluster_blast, taxa_colors = load_data()

## Data preprocessing

# Keep BGC that have a ClusterBlast similarity score with antismash DB greater than X%
cluster_blast_df = cluster_blast.copy()
cluster_blast_df['MAG'] = cluster_blast_df['sequence'].apply(lambda x : x.split("_")[0])
cluster_blast_df['sequence_w_type'] = cluster_blast_df['sequence'] + "_" + cluster_blast_df['cluster_type']


region_overview = region_summary.copy()
region_overview['MAG'] = region_overview['sequence'].apply(lambda x : x.split("_")[0])
region_overview['type'] = region_overview['type'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
region_overview = region_overview.explode('type').reset_index(drop=True)
region_overview = region_overview.drop_duplicates(subset='sequence')
region_overview = pd.merge(region_overview, virgo2_inventory[['MAG','FinalTaxonomy','classification']], on="MAG", how="left")
region_overview['Family'] = region_overview['classification'].apply(lambda x : x.split(";")[4][3:])
region_overview['Genus'] = region_overview['classification'].apply(lambda x : x.split(";")[5][3:])
region_overview['sequence_w_type'] = region_overview['sequence'] + "_" + region_overview['type']

# region_overview_filtered = region_overview[region_overview['sequence_w_type'].isin(cluster_blast_df[cluster_blast_df['similarity'] > 60]['sequence_w_type'].unique())].dropna()

# Make dictionary of colors
# colors = pc.qualitative.Set3 + pc.qualitative.Pastel1 + pc.qualitative.Set1 + pc.qualitative.Alphabet + pc.qualitative.Light24 + pc.qualitative.Prism + pc.qualitative.Antique + pc.qualitative.Pastel + pc.qualitative.Safe + pc.qualitative.Bold + pc.qualitative.Dark24 + pc.qualitative.Plotly + pc.qualitative.D3 + pc.qualitative.G10 + pc.qualitative.T10
custom_colors = (
    pc.qualitative.Plotly + pc.qualitative.Vivid + pc.qualitative.Set3 + pc.qualitative.Pastel1 + pc.qualitative.Set1 +
    pc.qualitative.Alphabet + pc.qualitative.Light24 + pc.qualitative.Prism + pc.qualitative.Antique + pc.qualitative.Pastel + 
    pc.qualitative.Safe + pc.qualitative.Bold + pc.qualitative.Dark24 + pc.qualitative.D3 + pc.qualitative.G10 +
    pc.qualitative.T10 + pc.qualitative.Dark2 + pc.qualitative.T10 +
    ["#E6194B", "#3CB44B", "#FFE119", "#0082C8", "#F58231","#911EB4", "#46F0F0", "#F032E6", "#D2F53C", "#008080","#AA6E28", "#FFFAC8", "#800000", "#AaffC3", "#808000","#FFD8B1", "#000080", "#808080", "#FFFFFF", "#000000","#B5651D", "#CFCFCF", "#4B0082", "#4682B4", "#D2691E","#FF69B4", "#CD5C5C", "#6A5ACD", "#708090", "#2E8B57"]
)
# custom_colors = pc.qualitative.Set3 + pc.qualitative.Set2 + pc.qualitative.Pastel1 + pc.qualitative.Light24
unique_cluster_type = region_overview['type'].unique()
color_mapping_type = {cluster: custom_colors[i % len(custom_colors)] for i, cluster in enumerate(unique_cluster_type)}
unique_cluster_clustertype = region_overview['most_similar_known_cluster_type'].unique()
color_mapping_clustertype = {cluster: custom_colors[i % len(custom_colors)] for i, cluster in enumerate(unique_cluster_clustertype)}
unique_cluster_compound = region_overview['most_similar_known_cluster'].unique()
color_mapping_compound = {cluster: custom_colors[i % len(custom_colors)] for i, cluster in enumerate(unique_cluster_compound)}


def display_barplot_bgc_taxonomic_level(annotation_column, top_value, threshold_similarity):
    # Function to generate grouped and sorted data
    def prepare_data(feature, top_value):
        df = region_overview.copy()
        if annotation_column == 'most_similar_known_cluster' or annotation_column == 'most_similar_known_cluster_type':
            df = df.dropna()
        df = df[df['sequence_w_type'].isin(cluster_blast_df[cluster_blast_df['similarity'] > threshold_similarity]['sequence_w_type'].unique())]
        df = df[[feature, annotation_column]]
        counts = df.groupby([feature, annotation_column]).size().unstack(fill_value=0)
        counts['Total'] = counts.sum(axis=1)
        counts = counts.sort_values('Total', ascending=False).head(top_value).drop(columns=['Total'])
        return counts

    genus_counts = prepare_data('Genus', top_value)

    if annotation_column == 'type':
        custom_colors = color_mapping_type
    elif annotation_column == 'most_similar_known_cluster_type':
        custom_colors = color_mapping_clustertype
    elif annotation_column == 'most_similar_known_cluster':
        custom_colors = color_mapping_compound

    # Make dictionary of colors
    legend_items = list(genus_counts.columns)
    legend_items = sorted(legend_items, key=str.lower)
    fig = go.Figure()
    for type_ in sorted(genus_counts.columns, key=str.lower):
        fig.add_trace(
            go.Bar(
                x=genus_counts[type_],
                y=genus_counts.index,
                name=type_,  # Legend item
                orientation="h",
                marker_color=custom_colors.get(type_, "#636efa"),
            )
        )

    # Update layout
    fig.update_layout(
        barmode="stack",
        title=f"{annotation_column} - Top {top_value}",
        xaxis_title="Count",
        yaxis_title="Genus",
        template="plotly_white",
        height=700,
        width=1400,  # Adjust width to fit single plot
    )

    fig.update_traces(marker=dict(line=dict(width=0)))
    st.plotly_chart(fig, use_container_width=True)

def get_all_taxa_region_table(taxa, feature, threshold = None):

    region = region_summary.copy()
    region['MAG'] = region['sequence'].apply(lambda x : x.split("_")[0])
    # explode 'type' column
    region['type'] = region['type'].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else x
    )
    region = region.explode('type').reset_index(drop=True)

    # if feature == 'type':
    region = region.drop_duplicates(subset='sequence')
    region['sequence_w_type'] = region['sequence'] + "_" + region['type']

    region = pd.merge(region, virgo2_inventory[['MAG','FinalTaxonomy','classification']], on="MAG", how="left")

    # Keep region that have a ClusterBlast similarity score with antismash DB greater than X%
    if threshold:
        region = region[region['sequence_w_type'].isin(cluster_blast_df[cluster_blast_df['similarity'] > threshold]['sequence_w_type'].unique())]
    
    region = region[region['FinalTaxonomy'].str.contains(taxa)]

    # Initialize an empty DataFrame to store results
    all_taxa_table = pd.DataFrame()

    # Loop through each unique taxa in `FinalTaxonomy`
    for taxa in region['FinalTaxonomy'].unique():
        # Filter the data for the current taxa
        region_taxa = region[region['FinalTaxonomy'] == taxa]
                
        # Calculate the table with total numbers and proportions
        region_taxa_table = (
            region_taxa[['FinalTaxonomy', feature]]
            .value_counts()
            .reset_index()
        )
        
        region_taxa_table = pd.merge(
            region_taxa_table,
            region_taxa[['MAG', feature]]
            .groupby(feature)
            .nunique()
            .rename(columns={"MAG": "Number_unique_MAG"})
            .reset_index(),
            on=feature,
            how='left'
        )

        region_taxa_table['N_mag_virgo2'] = virgo2_inventory[virgo2_inventory['FinalTaxonomy'] == taxa]['MAG'].nunique()
        
        # Calculate proportions for each cluster_type
        region_taxa_table['proportion_within_taxa'] = (region_taxa_table['Number_unique_MAG'] / region_taxa_table['N_mag_virgo2']).round(4)

        # Append the current taxa table to the overall table
        all_taxa_table = pd.concat([all_taxa_table, region_taxa_table])

    return all_taxa_table.reset_index(drop=True)


def display_barplot_per_species(df, title=None):

    feature = df.columns[1]
    if feature == 'type':
        color_mapping = color_mapping_type
        width = 1000
        height = 500
    elif feature == 'cluster_type':
        color_mapping = color_mapping_type
        width = 1000
        height = 500
    elif feature == 'most_similar_known_cluster_type':
        color_mapping = color_mapping_clustertype
        width = 1000
        height = 500
    elif feature == 'most_similar_known_cluster':
        color_mapping = color_mapping_compound
        width = 2000
        height = 500

    fig = px.bar(
        df.sort_values('FinalTaxonomy'),
        x="FinalTaxonomy",
        y="proportion_within_taxa",
        color=feature,
        # text='proportion_within_taxa',
        color_discrete_map=color_mapping,
        template='plotly_white',
        width=width,
        height=height,
        title=title,
    )

    # Prepare annotations for the top of each bar
    # annotations = []
    # final_taxonomy_groups = df.groupby("FinalTaxonomy")
    # for taxonomy, group in final_taxonomy_groups:
    #     # Get the total height of the stack for the current taxonomy
    #     total_height = group["proportion_within_taxa"].sum()
    #     # Get the `N_mag_virgo2` value for the taxonomy
    #     n_mag_value = group["N_mag_virgo2"].iloc[0] 
    #     annotations.append(
    #         dict(
    #             x=taxonomy,
    #             y=total_height,
    #             text=f"N: {n_mag_value}",
    #             showarrow=False,
    #             font=dict(size=10, color="black"),
    #             yanchor="bottom"
    #         )
    #     )

    # Update layout with annotations
    fig.update_layout(
        barmode="stack",
        # xaxis_tickangle=45,
        # annotations=annotations,
    )
    # fig.update_traces(marker=dict(line=dict(width=0)))  # Removes the border line around the bars

    # Show the figure
    st.plotly_chart(fig)


def scatter_w_barplot(column_label):
    type_to_color = taxa_colors.drop('Text', axis=1).set_index("Taxa")['Color'].to_dict()

    # Process region_summary file
    region_overview_mibig = region_summary.copy()
    region_overview_mibig['MAG'] = region_overview_mibig['sequence'].apply(lambda x: x.split("_")[0])
    # explode 'type' column
    region_overview_mibig['type'] = region_overview_mibig['type'].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else x
    )
    region_overview_mibig = region_overview_mibig.explode('type').reset_index(drop=True)
    region_overview_mibig = region_overview_mibig.drop_duplicates(subset='sequence')
    region_overview_mibig = pd.merge(region_overview_mibig, virgo2_inventory[['MAG','FinalTaxonomy','classification']], on="MAG", how="left")

    # Fill na
    region_overview_mibig = region_overview_mibig.dropna()

    # Create new columns
    bin_edges = np.arange(0, 110, 5)
    region_overview_mibig['similarity_bin'] = pd.cut(region_overview_mibig['similarity'], bins=bin_edges, right=False, labels=bin_edges[:-1])
    region_overview_mibig['Length'] = region_overview_mibig['From_To'].apply(lambda x: int(x.split("_")[1]) - int(x.split("_")[0]))

    available_colors = px.colors.qualitative.Plotly

    if column_label == 'type':
        custom_colors = color_mapping_type
    elif column_label == 'most_similar_known_cluster_type':
        custom_colors = color_mapping_clustertype
    elif column_label == 'most_similar_known_cluster':
        custom_colors = color_mapping_compound
    elif column_label == 'FinalTaxonomy':
        custom_colors = type_to_color

    def get_color(taxa, type_to_color, available_colors):
        if taxa in type_to_color:
            return type_to_color[taxa]  # Use assigned color
        else:
            return random.choice(available_colors)  # Assign a random color

    # Create a new color mapping for the `type` column
    region_overview_mibig['color'] = region_overview_mibig[column_label].map(lambda x: get_color(x, custom_colors, available_colors))

    # Data to plot
    region_overview_mibig_to_plot = region_overview_mibig.copy()
    stacked_data = region_overview_mibig_to_plot.groupby(['similarity_bin', column_label], observed=False).size().reset_index(name='count')

    # Create a subplot with two plots
    fig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.7, 0.3],
        specs=[[{"type": "scatter"}, {"type": "bar"}]]
    )

    # Add scatter plot traces by category for interactive legend
    for category in region_overview_mibig[column_label].unique():
        category_data = region_overview_mibig_to_plot[region_overview_mibig_to_plot[column_label] == category]
        fig.add_trace(
            go.Scatter(
                x=category_data['Length'],
                y=category_data['similarity'],
                mode='markers',
                name=category,
                marker=dict(color=category_data['color'].iloc[0]),
                legendgroup=category,
                showlegend=True  # Ensures legends are shown
            ),
            row=1, col=1
        )

    # Add stacked bar plot traces by category
    for category in stacked_data[column_label].unique():
        subset = stacked_data[stacked_data[column_label] == category]
        subset_color = region_overview_mibig[region_overview_mibig[column_label] == category]['color'].iloc[0]
        fig.add_trace(
            go.Bar(
                x=subset['count'],
                y=subset['similarity_bin'],
                name=category,
                orientation='h',
                marker=dict(color=subset_color, line=dict(width=0)),
                legendgroup=category,  # Links scatter and bar plot legends
                showlegend=False  # Avoid duplicate legends
            ),
            row=1, col=2
        )

    # Update axes and layout
    fig.update_xaxes(title_text="BGC Length", row=1, col=1)
    fig.update_xaxes(title_text="Count", row=1, col=2)
    fig.update_yaxes(title_text="Similarity (%)", range=[-5, 105], row=1, col=1)
    fig.update_yaxes(title_text="Similarity (%)", range=[-5, 105], row=1, col=2)
    fig.update_xaxes(
        # title_text="Similarity (%)",
        showline=True,
        linecolor="black",
        showgrid=True,
        gridcolor="lightgray",
        ticks="outside",
        tickcolor="black",
        tickwidth=2,
        row=1, col=1
    )

    fig.update_xaxes(
        # title_text="Similarity (%)",
        showline=True,
        linecolor="black",
        showgrid=True,
        gridcolor="lightgray",
        ticks="outside",
        tickcolor="black",
        tickwidth=2,
        row=1, col=2
    )
    fig.update_yaxes(
        title_text="Similarity (%)",
        showline=True,
        linecolor="black",
        showgrid=True,
        gridcolor="lightgray",
        ticks="outside",
        tickcolor="black",
        tickwidth=2,
        range=[-5, 105],
        row=1, col=1
    )

    fig.update_yaxes(
        title_text="Similarity (%)",
        showline=True,
        linecolor="black",
        showgrid=True,
        gridcolor="lightgray",
        ticks="outside",
        tickcolor="black",
        tickwidth=2,
        range=[-5, 105],
        row=1, col=2
    )
    fig.update_layout(
        barmode='stack',
        width=2000,
        height=600,
        showlegend=True,
        template='plotly_white'
    )

    # Show the plot
    st.plotly_chart(fig, use_container_width=True)


virgo2_family_genus = virgo2_inventory[['classification','FinalTaxonomy']].copy()
virgo2_family_genus['Family'] = virgo2_family_genus['classification'].apply(lambda x : x.split(";")[4][3:])
virgo2_family_genus['Genus'] = virgo2_family_genus['classification'].apply(lambda x : x.split(";")[5][3:])


def page():

    st.header("Table")
    st.dataframe(region_overview)

    st.header("BGC detection", divider = 'grey')
    feature_for_barplot = st.selectbox("Choose a feature", ("type", "most_similar_known_cluster_type", "most_similar_known_cluster"), key='taxonomic_level')
    # threshold_similarity = st.number_input("Threshold (cluster_blast similarity, %) ", min_value=0, max_value=100, value=0)
    col1, col2 = st.columns([4,1])
    with col1:
        display_barplot_bgc_taxonomic_level(feature_for_barplot, 15, threshold_similarity=0)
    with col2:
        st.subheader("Genus in VIRGO2")
        st.dataframe(virgo2_family_genus['Genus'].value_counts().reset_index())
    
    # st.header("MIBiG similarity score", divider = 'grey')
    # feature_w_taxa = st.radio("Choose a feature", ["type", "most_similar_known_cluster_type", "most_similar_known_cluster","FinalTaxonomy"], key='feature_w_taxa', index=0, horizontal=True)
    # scatter_w_barplot(feature_w_taxa)

    st.header("Species comparison", divider = 'grey')
    species = st.text_input("Grep...:")
    feature_for_species_barplot = st.radio("Choose a feature", ["type", "most_similar_known_cluster_type", "most_similar_known_cluster"], key='feature_for_species_barplot', index=0, horizontal=True)
    all_taxa_table = get_all_taxa_region_table(species, feature_for_species_barplot)
    display_barplot_per_species(all_taxa_table)

# Run the page function
if __name__ == "__main__":
    page()
