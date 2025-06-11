# Genesis: Protein Structure Prediction via Evolutionary Algorithm

**An Introduction to AI Term Project**

Genesis is a proof-of-concept project for my Introduction to AI course,
exploring an Evolutionary Algorithm (EA) approach to *ab initio* protein
structure prediction. It aims to create generalized estimators that can predict
a protein's 3D shape directly from its amino acid sequence.

**Why Genesis?**

Protein folding is incredibly complex, and finding the right 3D shape from a
sequence is computationally challenging. While modern tools like AlphaFold use
deep learning, they can be "black boxes" and need huge datasets. Genesis
explores heuristic methods as a transparent, flexible, and computationally
adaptable alternative, especially for novel proteins or when data is scarce. It
bridges classic optimization with modern structural biology, aiming for a
useful, understandable solution.

**How Genesis Works:**

Genesis dynamically generates a set of functions (f1, ..., fN). These functions
take a protein sequence and perform operations to give basic structural
predictions, acting like non-linear transformations. A "composite function" is
then built by combining these, and its "fitness" is calculated as the total
error (loss) over a database of known structures. Through evolution (mutating
the coefficients of these functions), a generalized estimator is created. This
estimator, once accurate enough, can then predict the 3D positions of a new
protein, given its precomputed functions.

**Running Genesis:**

1.  **Get the Code:**
Ensure you have the project files. The main code is in the `source` folder.

2.  **Set Up Your Environment:**
Open your terminal/command prompt.
* **Create a Python Virtual Environment:** `python -m venv venv`
* **Activate it:**
* Windows: `venv\Scripts\activate`
* macOS/Linux: `source venv/bin/activate`

3.  **Install Needed Libraries:**
Go into the `source` folder: `cd source`
Install dependencies: `pip install -r requirements.txt`

4.  **Run the Project:**
`python main.py`

**Important Note:** Genesis doesn't have a graphical user interface (GUI) yet.
To control what runs, you'll need to **open `main.py` and comment/uncomment
lines** as indicated by comments in the code. This allows you to selectively
run different parts of the protein prediction process.

**Challenges and Future:**

Protein folding is computationally demanding. Genesis, being a
proof-of-concept, has a lot of room for optimization, especially since some
calculations can grow exponentially (O(X^N)). This project focuses on
demonstrating the approach's feasibility and is a stepping stone for more
optimized solutions in the future.
