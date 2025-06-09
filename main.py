# main.py
import streamlit as st
import pandas as pd

st.set_page_config(layout='wide')

from pages_content import home, quality, taxa_comparison

# Set up sidebar navigation with "Home" as the default page
st.sidebar.title("Navigation")
page = st.sidebar.radio("Content", ["Home", "Quality", "Taxa comparison"], index=0, label_visibility='hidden')

st.sidebar.subheader("Contact")
st.sidebar.write("jholm@som.umaryland.edu")

# Display the selected page
if page == "Home":
    home.page()
elif page == "Quality":
    quality.page()
elif page == "Taxa comparison":
    taxa_comparison.page()  
