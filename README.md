# Human-Environment Stochastic Model

This repository contains two implementations of a human-environment interaction model for forest pest dynamics:

---

## 🧪 MATLAB Version

The original model and simulations (including Figure 2 from the paper) are implemented in MATLAB.

📂 Files:  
- `pest_model.m`  
- `*.m` helper scripts  
📄 Requirements: MATLAB 2015+  
📍 Location: Root folder

---

## 🧠 Python RL-style Version

The Python version rewrites the stochastic model in a modular, reinforcement-learning–friendly architecture using NumPy and Matplotlib.

📂 Folder: `python-code/`  
💻 Run with:

```bash
cd python-code
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python main.py
