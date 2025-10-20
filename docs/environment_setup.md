📄 docs/environment_setup.md 
# ⚙️ Environment Setup Guide This document describes how to install and activate the Conda environment for the **FluHub: Vaccine Design and Evolutionary Analysis Pipeline**. It ensures that all required tools and packages are installed with the correct versions for reproducibility. --- ## 🧩 1. Install Conda (if not already installed) If you don’t already have Conda on your system: ### 🪟 Windows Download and install **Miniconda** from: 👉 [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html) ### 🍎 macOS / 🐧 Linux Run in your terminal:
bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

Then restart your terminal and verify:

conda --version

🧪 2. Clone this repository
git clone https://github.com/<your-username>/Fluhub_Vaccine_design_and_evolutionary_analyses.git
cd Fluhub_Vaccine_design_and_evolutionary_analyses

🧱 3. Create the Conda environment

Run the following command from the repository root:

conda env create -f environment.yml


This will:

Create an isolated Conda environment named flu_design

Install all bioinformatics tools (MAFFT, IQ-TREE, CD-HIT, SeqKit)

Install all required Python packages with exact versions

🚀 4. Activate the environment

After successful installation:

conda activate flu_design


You should now see something like:

(flu_design) user@hostname:~/Fluhub_Vaccine_design_and_evolutionary_analyses$


This prefix (flu_design) means you’re working inside the correct environment.

🧠 5. Verify installation

Run the following commands to ensure all tools are available:

mafft --version
iqtree2 --version
cd-hit -h | head -n 1
seqkit version
python --version


You should see:

Tool	Example Output
MAFFT	v7.505
IQ-TREE	v2.2.6
CD-HIT	v4.7
SeqKit	v2.10.1
Python	3.11.13

Then test Python dependencies:

python -c "import Bio, pandas, numpy, matplotlib, seaborn, logomaker, scipy, tqdm; print('✅ Python packages loaded successfully!')"


If that prints the success message, your environment is fully working!

📦 6. (Optional) Export the exact environment snapshot

If you later add or upgrade packages, export the environment again to keep it up to date:

conda env export > environment.lock.yml


This captures the exact versions of everything you used — ideal for reproducibility or publication.

🧹 7. Deactivate and remove (if needed)

Deactivate:

conda deactivate
