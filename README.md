Chromatogram
-------------------
This pythoin script is used for parsing sanger sequencing results (abi file).
It's usefull to create chrmatogram when used for a paper.  

Example:  
![avatar](data/test.png "Example output")



Usage
============

In python console
```python
from chromatogram import Chromatogram
abi = Chromatogram("path/to/input.ab1", 
                    seq_range=(50, 100), # sequence range to show in the plot
                    rev_complement=False, # show reverse complement sequence and plot
                    figsize=(10,5)
                    )
abi.plot() # set filename="output.pdf" to save your figure
abi.plot(seq_range=(100,150)) # show sequence range from 100 to 150
```

Command line:
```bash
python chromatogram.py input.ab1 output.pdf seq_from seq_to
```

Dependency
=============

+ Python
   * Biopython
   * Matplotlib

