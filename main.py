# main.py
import streamlit as st
import pandas as pd

st.set_page_config(layout='wide')

from pages_content import home, quality, taxa_comparison

# Set up sidebar navigation with "Home" as the default page
st.sidebar.title("Navigation")
page = st.sidebar.radio("Content", ["Home", "BGC identification", "Taxanomic comparison"], index=0, label_visibility='hidden')
st.sidebar.divider()
st.sidebar.subheader("Contact")
st.sidebar.markdown("""[J B Holm Lab website](https://www.jbholmlab.org)""")
st.sidebar.write("jholm@som.umaryland.edu")

# Display the selected page
if page == "Home":
    home.page()
elif page == "BGC identification":
    quality.page()
elif page == "Taxanomic comparison":
    taxa_comparison.page()  
