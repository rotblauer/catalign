"""Streamlit dashboard application entry point."""

# This file is the entry point for the Streamlit dashboard
# Run with: streamlit run catalign/viz/dashboard_app.py
# Or use: catalign viz

from catalign.viz.dashboard import create_dashboard_content

if __name__ == "__main__":
    create_dashboard_content()
