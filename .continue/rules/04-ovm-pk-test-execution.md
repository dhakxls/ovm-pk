# OVM-PK Test Execution Rules

## Test Suite Execution Command

When asked to run tests, use these exact commands in sequence:

```bash
cd /home/martinvo/ovm-pk
python -m pytest -v tests/
```

## Expected Test Sequence

The test suite will execute in this specific order:
1. `test_fetchers.py` - Downloads protein/ligand files
2. `test_prep.py` - Prepares structures
3. `test_docking.py` - Runs docking
4. `test_pose_selection.py` - Selects best poses
5. `test_md_prepare.py` - Prepares MD system
6. `test_md_prod.py` - Runs production MD
7. `test_md_equil.py` - Runs equilibration
8. `test_md_analysis.py` - Analyzes results

## Output Locations

Results will be stored in:
- `data/output/run-{timestamp}/` - Main results
- `data/output/run-{timestamp}/logs/` - Log files
- `data/work/` - Intermediate files

## Test Configuration

Tests use settings from `configs/prod_test.yaml` including:
- Default protein: 5VCC
- Default ligand: ketoconazole
- Docking parameters
- MD simulation settings

## Direct Execution Command

When instructed "run tests" or similar, execute this exact command:

```python
import subprocess
import os

def run_ovm_pk_tests():
    os.chdir('/home/martinvo/ovm-pk')
    result = subprocess.run(['python', '-m', 'pytest', '-v', 'tests/'], 
                          capture_output=True, 
                          text=True)
    return result.stdout, result.stderr

# Execute tests
output, errors = run_ovm_pk_tests()
print(output)
if errors:
    print("Errors:", errors)
```

## Test Progress Monitoring

The test suite uses rich progress bars and logging. Monitor:
- Test execution order
- Progress bars for each stage
- Log output for detailed status
- Final results summary