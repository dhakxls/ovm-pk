# OVM-PK Product Requirements

## Core Functionality

### Physics-Aware Docking
- Replicate experimental affinity (Ki) and ΔG
- Support pluggable docking physics modules
- Enable head-to-head baseline comparisons

### Success Criteria
1. Replication accuracy:
   - |ΔG_pred − ΔG_exp| ≤ 1.5–2.0 kcal/mol median
   - Pearson r ≥ 0.5 between predicted and experimental ΔG

2. New physics improvements:
   - Statistically significant better correlation
   - Lower error vs baseline
   - Reproducible results

3. Operational excellence:
   - One-command end-to-end run
   - Rich logging and artifacts
   - Green test suite
   - Reproducible environments
   - Pluggable scorers/engines

## Technical Requirements

### Docking Physics
- Metal coordination (heme Fe–N/O)
- Polarization/induction effects
- Anisotropic VdW
- pH-dependent protonation
- Explicit/bridging waters
- Side-chain flexibility
- Microstate-aware scoring

### Energy Functions
```
E_total = w_vdw*E_LJ + w_elec*E_Coulomb + w_pol*E_pol + 
         w_coord*E_coord + w_desolv*E_GB + w_tors*E_lig_int
```

### Performance
- Milliseconds per pose evaluation
- Support for batch processing
- Optional GPU acceleration

## Implementation Requirements

### Code Organization
- Clear interface definitions
- Pluggable components
- Comprehensive test coverage
- Rich logging and metrics

### Configuration
- YAML-based configuration
- Reproducible runs
- Version control friendly

### Testing
- Unit tests for all components
- Integration tests for workflows
- Regression tests for physics
- Performance benchmarks

## Documentation Requirements

### API Documentation
- Clear interface specifications
- Usage examples
- Configuration schema
- Testing guidelines

### Scientific Documentation
- Method descriptions
- Physics equations
- Parameter explanations
- Validation metrics

### User Documentation
- Installation guide
- Usage tutorials
- Configuration examples
- Troubleshooting guide

## Deliverables

### Code
- Engine implementations
- Scorer implementations
- Benchmark runners
- Test suites

### Documentation
- API documentation
- Scientific methods
- User guides
- Test reports

### Artifacts
- Benchmark results
- Validation metrics
- Performance reports
- Configuration examples