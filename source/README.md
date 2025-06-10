# Genesis Project - Protein Structure Prediction using Evolutionary Algorithms

## Setup

1. Create a virtual environment:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Add your PDB files to the `data/` directory

4. Update the file paths in `main.py` to point to your PDB files

5. Run the program:
   ```bash
   python main.py
   ```

## Note

This is a demonstration of evolutionary algorithms applied to feature-to-coordinate mapping. 
For actual protein structure prediction, use specialized tools like AlphaFold or RoseTTAFold.
