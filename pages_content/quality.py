# pages_content/data_quality.py

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# # Load data
# virgo2_inventory = pd.read_csv("data/MAG_inventory_VIRGO2_021623_30Jul2024.txt", sep="\t")
# region_summary_original = pd.read_csv("data/csv_files/region_summary.csv")
# sequence_length = pd.read_csv("data/sequence_lengths.txt", sep="\t", header=None)

# Caching the data loading functions to speed up the Streamlit app
@st.cache_data
def load_data():
    # Load data
    virgo2_inventory = pd.read_csv("data/MAG_inventory_VIRGO2_021623_30Jul2024.txt", sep="\t")
    region_summary_original = pd.read_csv("data/region_summary.csv")
    sequence_length = pd.read_csv("data/sequence_lengths.txt.gz", sep="\t", header=None)
    
    return virgo2_inventory, region_summary_original, sequence_length

# Load the data using the cached function
virgo2_inventory, region_summary_original, sequence_length = load_data()

# Data processing
region_df = region_summary_original.copy()
region_df['MAG'] = region_df['sequence'].apply(lambda x: x.split("_")[0])

all_mags = virgo2_inventory['MAG'].unique()
mag_w_antismash_result = region_df['MAG'].unique()
mag_no_antismash_result = [i for i in all_mags if i not in mag_w_antismash_result]

antismash_status = pd.concat([
    pd.DataFrame({"MAG": mag_w_antismash_result, "status": 1}),
    pd.DataFrame({"MAG": mag_no_antismash_result, "status": 0})
], axis=0).sort_values("MAG", ascending=True)

sequences_data = sequence_length.copy().rename(columns={0: 'sequence', 1: 'length'})
sequences_data['MAG'] = sequences_data['sequence'].apply(lambda x: x.split("_")[0])
sequences_data = pd.merge(sequences_data, antismash_status, on="MAG", how="left")


# Merge and process data for display
stack_antismash_status = antismash_status.merge(virgo2_inventory[['MAG', 'FinalTaxonomy']], on='MAG', how='left')
status_counts = stack_antismash_status.groupby(['FinalTaxonomy', 'status']).size().unstack(fill_value=0)

# Sort FinalTaxonomy by total count in descending order
status_counts['Total'] = status_counts.sum(axis=1)
status_counts = status_counts.sort_values(by='Total', ascending=False).drop(columns='Total')
status_counts_long = status_counts.reset_index().melt(id_vars='FinalTaxonomy', var_name='status', value_name='count')


# Functions for displaying data
def display_antismash_status_pie(antismash_status):
    status_counts = antismash_status['status'].value_counts().reset_index()
    status_counts.columns = ['status', 'count']
    color_map = {1: "blue", 0: "red"}
    status_counts['color'] = status_counts['status'].map(color_map)
    
    fig = px.pie(
        status_counts, 
        values='count', 
        names='status', 
        color='status',
        title='Proportion of antiSMASH results (MAG)',
        color_discrete_map=color_map
    )
    st.plotly_chart(fig)

def display_numerical_feature_comparison(mag_inventory, antismash_status):
    numerical_columns = [col for col in mag_inventory.select_dtypes(include=['float64', 'int64']).columns 
                         if col not in ['Timepoint', 'FilterContam', 'red_value', 'warnings']]
    
    n_cols = 4
    n_rows = (len(numerical_columns) + n_cols - 1) // n_cols
    fig = make_subplots(rows=n_rows, cols=n_cols, subplot_titles=numerical_columns)

    merged_data = mag_inventory.merge(antismash_status, on="MAG", how="left")
    color_map = {0: "red", 1: "blue"}

    for i, col in enumerate(numerical_columns):
        row, col_pos = divmod(i, n_cols)
        for status in merged_data['status'].unique():
            filtered_data = merged_data[merged_data['status'] == status]
            fig.add_trace(
                go.Box(y=filtered_data[col], name=f"Status {status}", boxmean=True, marker_color=color_map[status]),
                row=row + 1, col=col_pos + 1
            )
        fig.update_yaxes(title_text="", row=row + 1, col=col_pos + 1)

    fig.update_layout(height=500 * n_rows, width=1500, showlegend=False)
    st.plotly_chart(fig)

def plot_mean_sequence_length(mean_length_data):

    color_map = {1: "blue", 0: "red"}
    fig = go.Figure()

    for status in mean_length_data['status'].unique():
        filtered_data = mean_length_data[mean_length_data['status'] == status]
        fig.add_trace(
            go.Histogram(
                x=filtered_data['length'],
                name=f"Status {status}",
                marker_color=color_map[status],
                nbinsx=200
            )
        )

    fig.update_layout(
        title="Mean Sequence Length per MAG",
        height=500, width=600,
        xaxis_title="Sequence Length",
        yaxis_title="Frequency",
        barmode='overlay',
        xaxis=dict(range=[-20000, 1000000]),
        showlegend=True
    )
    st.plotly_chart(fig)


def plot_number_of_sequences(count_length_data):

    color_map = {1: "blue", 0: "red"}
    fig = go.Figure()

    for status in count_length_data['status'].unique():
        filtered_data = count_length_data[count_length_data['status'] == status]
        fig.add_trace(
            go.Histogram(
                x=filtered_data['length'],
                name=f"Status {status}",
                marker_color=color_map[status],
                nbinsx=200
            )
        )

    fig.update_layout(
        title="Number of Sequences per MAG",
        height=500, width=600,
        xaxis_title="Number of Sequences",
        yaxis_title="Frequency",
        barmode='overlay',
        xaxis=dict(range=[-1000, 10000]),
        showlegend=True
    )
    st.plotly_chart(fig)


def display_taxa_processed(taxa_filter=None):
    if taxa_filter is not None :
        to_plot = status_counts_long[status_counts_long['FinalTaxonomy'].str.contains(taxa_filter, case=False, na=False)] 
    else :
        to_plot = status_counts_long

    color_map = {1: "blue", 0: "red"}
    fig = px.bar(
        to_plot, 
        x='count', 
        y='FinalTaxonomy', 
        color='status', 
        title="Stacked Bar Plot of Status Counts by Final Taxonomy (Sorted by Total Count)", 
        labels={"FinalTaxonomy": "Final Taxonomy", "count": "Count", "status": "Status"}, 
        height=1000,
        text_auto=True,
        color_discrete_map=color_map
    )
    # fig.update_layout(
    #     xaxis_title="Final Taxonomy", yaxis_title="Count", xaxis_tickangle=45, barmode='stack',
    #     legend_title="Status", template="plotly_dark",
        
    # )
    # fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)

    st.plotly_chart(fig)

# Streamlit page function
def page():
    st.title("Data Quality")

    st.subheader("VIRGO2 inventory", divider='grey')
    st.dataframe(pd.merge(virgo2_inventory, antismash_status, on='MAG', how='left'))

    st.subheader("Proportion of antiSMASH results", divider='grey')
    col1, col2 = st.columns(2)
    with col1:
        # Plot on the top
        display_antismash_status_pie(antismash_status)

        # Plot on the botton
        taxa_selection = st.selectbox("Taxa selection", sorted(stack_antismash_status['FinalTaxonomy'].unique()))
        filtered_data = stack_antismash_status[stack_antismash_status['FinalTaxonomy'] == taxa_selection]
        status_counts = filtered_data.groupby('status').size().reset_index(name='count')
        color_map = {1: "blue", 0: "red"}
        status_counts['color'] = status_counts['status'].map(color_map)

        fig = px.pie(status_counts, values='count', names='status', color='status', color_discrete_map=color_map)
        st.plotly_chart(fig)

    with col2:
        taxa_filter = st.text_input("Grep...")
        display_taxa_processed(taxa_filter)

    zero_only = []
    one_only = []
    for i in status_counts_long['FinalTaxonomy'].unique():
        df_temp = status_counts_long[status_counts_long['FinalTaxonomy'] == i]
        if df_temp.loc[df_temp['count'] == 0].shape[0] != 0 :
            zero_count = df_temp.loc[df_temp['count'] == 0]
            if zero_count['status'].values[0] == 0:
                one_only.append(i)
            else :
                zero_only.append(i)

    # col1, col2 = st.columns(2)
    # with col1:
    #     st.subheader("No results at all", divider='grey')
    #     df_zero_only = virgo2_inventory[virgo2_inventory['FinalTaxonomy'].isin(zero_only)]['FinalTaxonomy'].value_counts().reset_index().sort_values("FinalTaxonomy")
    #     # st.dataframe(pd.DataFrame(sorted(zero_only)))
    #     st.dataframe(df_zero_only.reset_index(drop=True))

    # with col2:
    #     st.subheader(f"100% of results", divider='grey')
    #     df_one_only = virgo2_inventory[virgo2_inventory['FinalTaxonomy'].isin(one_only)]['FinalTaxonomy'].value_counts().reset_index().sort_values("FinalTaxonomy")
    #     st.dataframe(df_one_only.reset_index(drop=True))
    
    # # with col3:
    #     st.subheader("All taxa", divider='grey')
    #     st.dataframe(status_counts_long.sort_values(by="FinalTaxonomy",ascending=True).reset_index(drop=True))

    st.subheader("Boxplot of numerical metadata", divider='grey')
    st.dataframe(pd.DataFrame(virgo2_inventory.isna().sum()[virgo2_inventory.isna().sum() != 0], columns=['NaN']).transpose())
    display_numerical_feature_comparison(virgo2_inventory, antismash_status)

    st.subheader("Distribution of sequences lengths", divider='grey')
    mean_length_data = sequences_data[['length', 'MAG']].groupby("MAG").mean().reset_index().merge(antismash_status, on="MAG")
    count_length_data = sequences_data[['length', 'MAG']].groupby("MAG").count().reset_index().merge(antismash_status, on="MAG")
    col1, col2 = st.columns(2)
    with col1:
        plot_mean_sequence_length(mean_length_data)
    with col2:
        plot_number_of_sequences(count_length_data)

# Run the page function
if __name__ == "__main__":
    page()