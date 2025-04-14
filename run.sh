#!/bin/bash 
python -m streamlit run streamlit_app.py \ --server.port 8000 \ --server.address 0.0.0.0 \ --server.enableCORS false \ --server.enableXsrfProtection false \ --server.enableWebsocketCompression false \ --server.enableWebsocketConnectionReuse false \ --server.headless true \ --browser.serverAddress "https://gist-b.azurewebsites.net" \ --server.runOnSave false \ --client.showErrorDetails false \ --client.toolbarMode minimal \ --server.enableStaticServing true


