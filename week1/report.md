# Deliverable 1 Report  

Overall the Python code was converted sucessfully and does run in Codon. However some outputs vary depending on whether the program is ran with Python or Codon. For data1 and data2 tests both Codon and Python have the same outputs however for data3 and data4 tests show different outputs depending on whether it is ran in Python or Codon. For N50 calculations we have data1-3 outputs being the same when ran in Codon and Python but N50 for data4 when ran in Codon differs to what is output when ran in Python. After alot of attempts to fix this it seems that this may be due to Codon’s built-in sets/lists not guaranteeing the same iteration order as CPython, which would be affecting tie-breaking in graph traversal. Ended up partially solving this by forcing CPython set semantics via interop wrappers, but this only made data2 tests 100% the same across Python and Codon. Data3 and data4 output saw improvements with this fix but still was not 100% the same across Python and Codon versions. Any advise on how to develope a proper work around for this issue would be appreciated.  

Will walk throught the general steps taken below.  

## Steps:  

1. Cloning and setup  
   - Cloned the genome-assembly repository and explored the code.  
   - Verified that the provided Python version (main.py, dbg.py, utils.py) ran correctly on the four datasets data1 to data4.  
  
2. Baseline run in Python  
   - Ran python3 main.py ../data/dataX for each dataset.  
   - was unable to recreeate the full table's values from README/report in genome-assembly repo as not all tools/info to do was available.  

3. Codon conversion  
   - Converted the genome assembler to run under Codon.  
   - Replaced Python-specific parts (e.g., os.path.join) with string concatenation that worked in Codon since os.path.join was not working.  
   - Added Codon interop for data structures (e.g., CPython set helpers) so that iteration order and hashing matched CPython semantics as best they could which lead to an improvement in all the data1-4's outputs. (still isn't perfect)  
   - Ensured dbg.py was Codon compatible while preserving the same algorithm as the original Python version.  

4. Automation  
   - Wrote tests.sh to run both Python and Codon versions across all datasets and compare outputs.  
   - Wrote evaluate.sh to benchmark runtime and compute N50 values automatically and it produces a formatted results table.  
        - Time format for the table is hours : minutes : seconds.  
  
5. CI/CD  
   - Finally, added GitHub Actions workflow (.github/workflows/ci.yml) to automatically install Codon, set CODON_PYTHON, and run the evaluation script on every push.  

## "Gotchas"  

- Matplotlib: Original code imported matplotlib but did not use it in the final version. This had to be handled carefully since Codon cannot run matplotlib.  
- Recursion depth: Dataset 4 data4 required increasing the stack size ulimit -s 8192000 to avoid recursion limit errors for deep graphs.  
- Codon type inference: Codon sometimes mis-inferred types for empty dicts and sets, leading to errors like 'NoneType' object has no attribute '__hash__'. My best fix for this was by using CPython interop @python helpers so Codon defers to real Python sets/dicts.  
- Iteration order: After looking into the issue for a while, and to the best of my understanding it seems Codon’s built-in sets/lists did not guarantee the same iteration order as CPython, which affected tie-breaking in graph traversal. Partially solved this by forcing CPython set semantics via interop wrappers.  
- Outputs: Some Codon outputs initially differed from Python outputs (especially in datasets 2–4). Through debugging, wit was found it was due to type inference and ordering issues. After fixing with interop sets, Codon and Python results aligned for data1 and data2 with improvements but still flaws in data3 and data4 outputs.  

To run evaluate.sh or tests.sh please run these commands in this order:  
   chmod +x evaluate.sh  
   ./evaluate.sh  

   or  

   chmod +x tests.sh  
   ./tests.sh  

To run the code manually without evaluation.sh or without tests.sh please be sure to run the following commands in the order shown below:
   python3 main.py ../data/data1  
   python3 main.py ../data/data2  
   python3 main.py ../data/data3  
   ulimit -s 8192000  
   python3 main.py ../data/data3  
   
   export CODON_PYTHON=/usr/lib/x86_64-linux-gnu/libpython3.12.so.1.0  
   codon run -release main.py ../data/data1  
   codon run -release main.py ../data/data2  
   codon run -release main.py ../data/data3  
   ulimit -s 8192000  
   codon run -release main.py ../data/data4  
 


