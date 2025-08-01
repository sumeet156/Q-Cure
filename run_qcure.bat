@echo off
echo Starting Q-Cure Quantum Drug Interaction Predictor...
echo.
cd /d "%~dp0\app"
streamlit run streamlit_app.py --server.address localhost --server.port 8501
pause
